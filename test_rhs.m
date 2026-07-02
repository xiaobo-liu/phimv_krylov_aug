function result = test_rhs(opts)
%TEST_RHS Run right-hand-side beta/delta sensitivity experiments.

if nargin < 1 || ~isstruct(opts)
    error('test_rhs:MissingOptions', ...
        'test_rhs requires an options struct from run_all_test.');
end

required_fields = {'matnames', 'fig_dir', 'data_dir', 'n', 's', ...
    'm_max', 'matrix_scale', 'fixed_beta', 'fixed_delta', ...
    'delta_values', 'beta_values', 'rerun_test'};
missing = required_fields(~isfield(opts, required_fields));
if ~isempty(missing)
    error('test_rhs:MissingOptions', ...
        'test_rhs options must be specified in run_all_test. Missing: %s.', ...
        strjoin(missing, ', '));
end

%% test parameters from run_all_test

matnames = opts.matnames;
num_matnames = length(matnames);

fig_dir = opts.fig_dir;
data_dir = opts.data_dir;
rerun_test = opts.rerun_test;

n = opts.n;
s = opts.s;
m_max = opts.m_max;
matrix_scale = opts.matrix_scale;
fixed_beta = opts.fixed_beta;
fixed_delta = opts.fixed_delta;
delta_values = opts.delta_values;
beta_values = opts.beta_values;

%% run tests

main_loop = tic;

fprintf('Running right-hand-side sensitivity experiments\n');

% The output structure result contains one field per test matrix:
%   result.matrices.(matname)  delta and beta sensitivity data
%   result.elapsed_minutes     time spent by the sensitivity driver
result = struct();
result.matrices = struct();

for imat = 1:num_matnames
    matname = matnames{imat};

    [A, ~] = testmat(matname, n, matrix_scale);
    fprintf('Running test for matrix: %s\n', matname);
    % Use one family per matrix so regenerated cases change only beta or delta.
    rhs_family = make_rhs_family(size(A, 1), s);

    fprintf('Delta sensitivity test for matrix: %s\n', matname);
    delta_results = cell(1, numel(delta_values));
    for k = 1:numel(delta_values)
        beta = fixed_beta;
        delta = delta_values(k);

        dataname = fullfile(data_dir, sprintf( ...
            'result_delta_%s_beta_%.4g_delta_%.4g.mat', ...
            matname, beta, delta));
        case_result = load_or_run_case(dataname, rerun_test, ...
            matname, n, s, m_max, beta, delta, matrix_scale, rhs_family);

        case_result.sens_name = 'delta';
        case_result.sens_value = case_result.delta;
        delta_results{k} = case_result;

        fprintf(['  beta=%8.2e  delta=%8.2e  ||B||_2=%8.2e  ', ...
            'cond(X_K)=%8.2e  sigma(X_K)=[%8.2e,%8.2e]\n'], ...
            case_result.beta, case_result.delta, ...
            case_result.norm_B, case_result.cond_XK, ...
            case_result.sigma_min_XK, case_result.sigma_max_XK);
    end

    fig = plot_fov_sens(delta_results, 'delta');
    fov_delta_file = fullfile(fig_dir, sprintf('%s_fov_delta_sens', matname));
    save_test_fig(fig, fov_delta_file);
    % Leave figures open so run_all_test shows the generated plots.
    % close(fig)

    fig = plot_bound_sens(delta_results, 'delta');
    conv_delta_file = fullfile(fig_dir, ...
        sprintf('%s_convergence_delta_sens', matname));
    save_test_fig(fig, conv_delta_file);
     % close(fig)

    fprintf('Beta sensitivity test for matrix: %s\n', matname);
    beta_results = cell(1, numel(beta_values));
    for k = 1:numel(beta_values)
        beta = beta_values(k);
        delta = fixed_delta;

        dataname = fullfile(data_dir, sprintf( ...
            'result_beta_%s_beta_%.4g_delta_%.4g.mat', ...
            matname, beta, delta));
        case_result = load_or_run_case(dataname, rerun_test, ...
            matname, n, s, m_max, beta, delta, matrix_scale, rhs_family);

        case_result.sens_name = 'beta';
        case_result.sens_value = case_result.beta;
        beta_results{k} = case_result;

        fprintf(['  beta=%8.2e  delta=%8.2e  ||B||_2=%8.2e  ', ...
            'cond(X_K)=%8.2e  sigma(X_K)=[%8.2e,%8.2e]\n'], ...
            case_result.beta, case_result.delta, ...
            case_result.norm_B, case_result.cond_XK, ...
            case_result.sigma_min_XK, case_result.sigma_max_XK);
    end

    fig = plot_fov_sens(beta_results, 'beta');
    fov_beta_file = fullfile(fig_dir, sprintf('%s_fov_beta_sens', matname));
    save_test_fig(fig, fov_beta_file);
    % close(fig)

    fig = plot_bound_sens(beta_results, 'beta');
    conv_beta_file = fullfile(fig_dir, ...
        sprintf('%s_convergence_beta_sens', matname));
    save_test_fig(fig, conv_beta_file);
    % close(fig)
    
    % Store the two one-parameter sweeps separately for this matrix.
    matrix_result = struct();
    matrix_result.effective_n = size(A, 1);
    matrix_result.delta = [delta_results{:}];
    matrix_result.beta = [beta_results{:}];
    result.matrices.(matname) = matrix_result;
end

result.elapsed_minutes = toc(main_loop) / 60;

fprintf('Producing the right-hand-side sensitivity results took %.2f minutes.\n', ...
    result.elapsed_minutes);

end

function case_result = load_or_run_case(dataname, rerun_test, ...
    matname, n, s, m_max, beta, delta, matrix_scale, rhs_family)
%LOAD_OR_RUN_CASE Load stored data, regenerating when requested or missing.

if rerun_test || ~exist(dataname, 'file')
    case_result = fov_krylov_comparison(matname, n, s, m_max, ...
        beta, delta, matrix_scale, rhs_family);
    result_file = struct('result', case_result);
    save(dataname, '-struct', 'result_file');
else
    result_file = load(dataname, 'result');
    case_result = result_file.result;
end

end
