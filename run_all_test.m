% Run all field-of-values and Krylov convergence tests.

warning off
rng default
format compact

code_dir = fileparts(mfilename('fullpath'));
if ~isempty(code_dir), cd(code_dir); end
run_tests = struct();

% If false, current stored data are loaded and the plots are rebuilt.
run_tests.examples = false; % "run_tests.examples = true" takes about 30 minutes
run_tests.rhs = false;      % "run_tests.rhs = true" takes about 120 minutes

n = 40;
m_max = 30;
s = 5;

%% Parameters

shared_opts = struct();
shared_opts.data_dir = fullfile(code_dir, 'data');
shared_opts.fig_dir = fullfile(code_dir, 'fig');
shared_opts.n = n;
shared_opts.s = s;                   % compute phi_0,...,phi_s
shared_opts.m_max = m_max;           % maximal Krylov dimension
shared_opts.matrix_scale = 6;        % keeps expm(K) moderate

% Five examples
example_opts = shared_opts;
example_opts.mats = {'poisson', 'kms', 'grcar', 'dorr', 'triw'};
example_opts.rerun_test = run_tests.examples;
example_opts.beta = 10;               % common right-hand-side scale
example_opts.delta = 0.1;             % perturbation level in q + delta*r_j

% Right-hand-side sensitivity tests
rhs_opts = shared_opts;
rhs_opts.matnames = {'poisson', 'grcar'};
rhs_opts.fixed_beta = 10;
rhs_opts.fixed_delta = 0.1;
% Vary one right-hand-side parameter at a time.
rhs_opts.delta_values = [0, 1e-3, 1e-1, 10];
rhs_opts.beta_values = [1, 4, 10, 25];
rhs_opts.rerun_test = run_tests.rhs;

% output directories
output_dirs = {shared_opts.fig_dir, shared_opts.data_dir};
for k = 1:numel(output_dirs)
    if ~exist(output_dirs{k}, 'dir'), mkdir(output_dirs{k}); end
end

%% run selected tests

main_loop = tic;

fprintf('Running all tests\n');
fprintf('  n             = %d\n', shared_opts.n);
fprintf('  s             = %d\n', shared_opts.s);
fprintf('  m_max         = %d\n', shared_opts.m_max);
fprintf('  matrix_scale  = %.2e\n', shared_opts.matrix_scale);
fprintf('  data_dir      = %s\n', shared_opts.data_dir);
if example_opts.rerun_test
    example_mode = 'regenerate data and plots';
else
    example_mode = 'load stored data and rebuild plots';
end
if rhs_opts.rerun_test
    rhs_mode = 'regenerate data and plots';
else
    rhs_mode = 'load stored data and rebuild plots';
end
fprintf('Parameters in five examples (%s)\n', example_mode);
fprintf('  example_mats  = %s\n', strjoin(example_opts.mats, ', '));
fprintf('  beta          = %.2e\n', example_opts.beta);
fprintf('  delta         = %.2e\n', example_opts.delta);

% Convert the numeric arrays into strings for screen printing.
delta_parts = arrayfun(@(v) sprintf('%.4g', v), ...
    rhs_opts.delta_values, 'UniformOutput', false);
beta_parts = arrayfun(@(v) sprintf('%.4g', v), ...
    rhs_opts.beta_values, 'UniformOutput', false);
fprintf('Parameters in right-hand-side sensitivity tests (%s)\n', rhs_mode);
fprintf('  rhs_matrices  = %s\n', strjoin(rhs_opts.matnames, ', '));
fprintf('  fixed beta    = %.2e\n', rhs_opts.fixed_beta);
fprintf('  delta values  = [%s]\n', strjoin(delta_parts, ', '));
fprintf('  fixed delta   = %.2e\n', rhs_opts.fixed_delta);
fprintf('  beta values   = [%s]\n', strjoin(beta_parts, ', '));

% The output structure result records the run setup and all loaded or
% regenerated data:
%   result.settings.shared     common dimensions and output directories
%   result.settings.examples   parameters for the five model examples
%   result.settings.rhs        parameters for the sensitivity experiments
%   result.examples            five example results and timing
%   result.rhs                 right-hand-side sensitivity results and timing
%   result.elapsed_minutes     total run time
result = struct();
result.settings = struct();
result.settings.run_tests = run_tests;
result.settings.shared = shared_opts;

result.settings.examples = struct();
result.settings.examples.mats = example_opts.mats;
result.settings.examples.rerun_test = example_opts.rerun_test;
result.settings.examples.beta = example_opts.beta;
result.settings.examples.delta = example_opts.delta;

result.settings.rhs = struct();
result.settings.rhs.matnames = rhs_opts.matnames;
result.settings.rhs.fixed_beta = rhs_opts.fixed_beta;
result.settings.rhs.fixed_delta = rhs_opts.fixed_delta;
result.settings.rhs.delta_values = rhs_opts.delta_values;
result.settings.rhs.beta_values = rhs_opts.beta_values;
result.settings.rhs.rerun_test = rhs_opts.rerun_test;

result.examples = test_examples(example_opts);
result.rhs = test_rhs(rhs_opts);

result.elapsed_minutes = toc(main_loop) / 60;

fprintf('All tests took %.2f minutes.\n', result.elapsed_minutes);
fprintf('Results stored in result.examples and result.rhs.\n');
