% Test field-of-values estimates and Krylov convergence

code_dir = fileparts(mfilename('fullpath'));
if ~isempty(code_dir), cd(code_dir); end
root_dir = fileparts(code_dir);
warning off

rng default
format compact

fig_dir = fullfile(root_dir, 'fig');
data_dir = fullfile(code_dir, 'data');
dirs = {fig_dir, data_dir};
for k = 1:numel(dirs)
    if ~exist(dirs{k}, 'dir'), mkdir(dirs{k}); end
end

%% test parameters

mats = {'laplace', 'kms', 'grcar', 'dorr', 'triw'};
num_mats = length(mats);

% rerun the test or load from existing test data
rerun_test = false;

n = 40;
s = 5;                   % compute phi_0,...,phi_s
m_max = 30;              % maximal Krylov dimension

rhs_scale = 10;          % large values enlarge the Euclidean F(K) estimate
noise_scale = 1e-2;        % small perturbation of the common RHS direction
matrix_scale = 6;        % keeps expm of the augmented matrices moderate

%% run all tests

main_loop = tic;
results = struct();

for k = 1:num_mats
    matname = mats{k};
    fprintf('Running test for matrix: %s\n', matname);
    
    dataname = fullfile(data_dir, sprintf('result_%s.mat', matname));
    if rerun_test
        result = run_mat_comparison(matname, n, s, m_max, rhs_scale, noise_scale, matrix_scale);
        results.(matname) = result;
        save(dataname, 'result');
    else
        load(dataname, 'result');
    end
    
    figs = plot_result(result);

    save_test_fig(figs.fov, fullfile(fig_dir, sprintf('comparison_fov_%s', matname)));
    save_test_fig(figs.error, fullfile(fig_dir, sprintf('comparison_error_%s', matname)));

    fprintf('  ||W Q_K - Q_K K|| / ||W Q_K||   = %.2e\n', result.res_K);
    fprintf('  ||W Q_X - Q_X K_X|| / ||W Q_X|| = %.2e\n', result.res_X);
    fprintf('  max |err_W - err_X|             = %.2e\n', result.err_gap);
    fprintf('  ||B||_2                         = %.2e\n', result.norm_B);
    fprintf('  cond(X_K)                       = %.2e\n', result.cond_XK);
end

fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
fprintf('Figures were saved in %s.\n', fig_dir);
