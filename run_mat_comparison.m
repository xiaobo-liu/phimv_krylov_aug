function result = run_mat_comparison(matname, n, s, m_max, rhs_scale, noise_scale, matrix_scale)
%RUN_MAT_COMPARISON Run one FOV/Krylov comparison for one test matrix.

num_theta = 200;         % number of support lines for FOV plots

[A, mat_label] = testmat(matname, n, matrix_scale);
n = size(A, 1);
U = rhs_vectors(n, s, rhs_scale, noise_scale); % U = [b0,b1,...,bs].

% Build K, W, the KIOPS basis X_K, and the orthonormal K_X system.
aug = build_aug_mat(A, U);

max_krylov_dim = size(aug.K,1) - 1;
m_max = min(m_max, max_krylov_dim);
mvals = 1:m_max;

% Reference projected vector y = Pi_K exp(K)c for all three Krylov methods 
z_ref = expm(full(aug.K)) * aug.c;
y_ref = z_ref(1:n);
norm_y = max(norm(y_ref), realmin);

err_kiops = zeros(size(mvals));
err_W = zeros(size(mvals));
err_X = zeros(size(mvals));
for j = 1:length(mvals)
    m = mvals(j);
    [y, ~] = phimv_kry_aug(A, U, s, m, 0, 'kiops');
    err_kiops(j) = norm(y-y_ref) / norm_y;

    [y, ~] = phimv_kry_aug(A, U, s, m, 0, 'block');
    err_W(j) = norm(y-y_ref) / norm_y;

    [y, ~] = phimv_kry_aug(A, U, s, m, 0, 'orth');
    err_X(j) = norm(y-y_ref) / norm_y;
end

range_opt.num_theta = num_theta;
range_K = fov_boundary(aug.K, range_opt);
range_MK = fov_boundary(aug.RK_K_RKinv, range_opt);
range_W = fov_boundary(aug.W, range_opt);
range_A = fov_boundary(A, range_opt);
range_J = fov_boundary(aug.Js, range_opt);

% The K bound contains ||B||/2, while W contains 1/2.
bnd_K = fov_boundary_bnd([range_A, range_J], 0.5*norm(aug.B), num_theta);
bnd_W = fov_boundary_bnd([range_A, range_J], 0.5, num_theta);

crouzeix = 1 + sqrt(2);
est_K  = converg_bestpoly_bnd(range_K,  norm(aug.c), norm_y, mvals, crouzeix);
est_MK = converg_bestpoly_bnd(range_MK, norm(aug.b), norm_y, mvals, crouzeix);
est_W  = converg_bestpoly_bnd(range_W,  norm(aug.b), norm_y, mvals, crouzeix);

% Residual test for the identities W*Q = Q*K.
res_K = norm(full(aug.W*aug.QK - aug.QK*aug.K), 'fro') / ...
    max(norm(full(aug.W*aug.QK), 'fro'), realmin);
res_X = norm(full(aug.W*aug.QX - aug.QX*aug.KX), 'fro') / ...
    max(norm(full(aug.W*aug.QX), 'fro'), realmin);

result.matname = matname;
result.mat_label = mat_label;
result.n = n;
result.s = s;
result.mvals = mvals;

result.err_kiops = err_kiops;
result.err_W = err_W;
result.err_X = err_X;
result.est_K = est_K;
result.est_MK = est_MK;
result.est_W = est_W;

result.range_K = range_K;
result.range_MK = range_MK;
result.range_W = range_W;
result.range_A = range_A;
result.bnd_K = bnd_K;
result.bnd_W = bnd_W;

result.eig_K = eig(full(aug.K));

result.res_K = res_K;
result.res_X = res_X;
result.err_gap = max(abs(err_W-err_X));
result.cond_XK = cond(full(aug.XK));
result.norm_B = norm(aug.B);

end

function U = rhs_vectors(n, s, rhs_scale, noise_scale)
%RHS_VECTORS  Right-hand side vectors U = [b0,b1,...,bs].

U = zeros(n, s+1);

b0 = randn(n, 1);
U(:,1) = b0 / norm(b0);

base = randn(n, 1);
base = base / norm(base);

for j = 1:s
    noise = randn(n, 1);
    noise = noise / norm(noise);

    % Common direction with a small independent perturbation.
    bj = base + noise_scale*noise;
    U(:,j+1) = rhs_scale * bj / norm(bj);
end

end
