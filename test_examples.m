function result = test_examples(opts)
%TEST_EXAMPLES Run the five model matrix examples.

if nargin < 1 || ~isstruct(opts)
    error('test_examples:MissingOptions', ...
        'test_examples requires an options struct from run_all_test.');
end

required_fields = {'mats', 'fig_dir', 'data_dir', 'rerun_test', ...
    'n', 's', 'm_max', 'beta', 'delta', 'matrix_scale'};
missing = required_fields(~isfield(opts, required_fields));
if ~isempty(missing)
    error('test_examples:MissingOptions', ...
        'test_examples options must be specified in run_all_test. Missing: %s.', ...
        strjoin(missing, ', '));
end

%% test parameters from run_all_test

mats = opts.mats;
num_mats = length(mats);

fig_dir = opts.fig_dir;
data_dir = opts.data_dir;
rerun_test = opts.rerun_test;

n = opts.n;
s = opts.s;
m_max = opts.m_max;
beta = opts.beta;
delta = opts.delta;
matrix_scale = opts.matrix_scale;

%% run tests

main_loop = tic;

fprintf('Running five-example experiments\n');

% The output structure result contains one field per model matrix:
%   result.cases.(matname)  loaded or regenerated comparison data
%   result.elapsed_minutes  time spent by the five-example driver
result = struct();
result.cases = struct();

for k = 1:num_mats
    matname = mats{k};
    fprintf('Running test for matrix: %s\n', matname);

    dataname = fullfile(data_dir, sprintf('result_%s.mat', matname));
    case_result = load_or_run_case(dataname, rerun_test, ...
        matname, n, s, m_max, beta, delta, matrix_scale);

    figs = plot_matrix(case_result);
    save_test_fig(figs.fov, ...
        fullfile(fig_dir, sprintf('comparison_fov_%s', matname)));
    save_test_fig(figs.error, ...
        fullfile(fig_dir, sprintf('comparison_error_%s', matname)));
    if isfield(figs,'error_sym')
        save_test_fig(figs.error_sym, ...
            fullfile(fig_dir, sprintf('comparison_error_%s_sym', matname)));
    end
    % Leave figures open so run_all_test shows the generated plots.
    % close(figs.fov)
    % close(figs.error)

    result.cases.(matname) = case_result;

    fprintf('  ||W Q_K - Q_K K|| / ||W Q_K||   = %.2e\n', case_result.res_K);
    fprintf('  ||W Q_X - Q_X K_X|| / ||W Q_X|| = %.2e\n', case_result.res_X);
    fprintf('  max |err_W - err_X|             = %.2e\n', case_result.err_gap);
    fprintf('  ||B||_2                         = %.2e\n', case_result.norm_B);
    fprintf('  cond(X_K)                       = %.2e\n', case_result.cond_XK);
end

result.elapsed_minutes = toc(main_loop) / 60;

fprintf('Producing the example results took %.2f minutes.\n', ...
    result.elapsed_minutes);

end

function case_result = load_or_run_case(dataname, rerun_test, ...
    matname, n, s, m_max, beta, delta, matrix_scale)
%LOAD_OR_RUN_CASE Load stored data, regenerating when requested or missing.

if rerun_test || ~exist(dataname, 'file')
    case_result = fov_krylov_comparison(matname, n, s, m_max, ...
        beta, delta, matrix_scale);
    result_file = struct('result', case_result);
    save(dataname, '-struct', 'result_file');
else
    result_file = load(dataname, 'result');
    case_result = result_file.result;
end

end
