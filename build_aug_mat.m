function aug = build_aug_mat(A, U, basis)
%BUILD_AUG_MAT Build augmented matrices.
%
%   aug = build_aug_mat(A,U) builds all matrices needed by the
%   numerical experiment.
%
%   aug = build_aug_mat(A,U,basis) builds only the requested realization.
%   Accepted basis names are 'kiops'/'k', 'block'/'w', or 'orth'/'kx'.

if nargin < 3 || isempty(basis)
    basis = 'all';
end
basis = lower(char(basis));

[n, cols] = size(U);
s = cols - 1;

b0 = U(:,1);
Js = spdiags(ones(s,1), 1, s, s);
aug.Js = Js;

build_all = strcmp(basis, 'all');
is_kiops_basis = any(strcmp(basis, {'kiops', 'k'}));
is_block_basis = any(strcmp(basis, {'block', 'w'}));
is_orth_basis = any(strcmp(basis, {'orth', 'kx'}));

build_kiops = build_all || is_kiops_basis;
build_orth = build_all || is_orth_basis;
build_block = build_all || is_block_basis || build_orth;

if ~(build_all || is_kiops_basis || is_block_basis || is_orth_basis)
    error(['Unknown basis "%s". Use "kiops"/"k", "block"/"w", ', ...
        '"orth"/"kx", or "all".'], basis);
end

if build_kiops
    % K is the KIOPS coordinate matrix.
    B = U(:,end:-1:2);
    es = zeros(s,1);
    es(s) = 1;
    aug.K = [A, B; sparse(s,n), Js];
    aug.c = [b0; es];
    aug.B = B;
end

if build_block
    % W is the large block matrix with fixed coupling E.
    bb = reshape(U(:,2:end), n*s, 1);
    J = kron(Js, speye(n));
    E = [speye(n), sparse(n, n*(s-1))];

    aug.W = [A, E; sparse(n*s,n), J];
    aug.b = [b0; bb];
    aug.J = J;
    aug.E = E;
    aug.bb = bb;
end

if build_orth
    % X_K = [J^(s-1)bb, ..., Jbb, bb] is the KIOPS lower basis.
    XK = zeros(n*s, s);
    for j = 1:s
        XK(:,j) = aug.J^(s-j) * aug.bb;
    end

    % Cholesky QR for orthonormal basis X for range(X_K).
    R = chol(XK'*XK);
    X = XK / R;

    % K_X is the reduced matrix for the orthonormal basis X.
    G = X' * (aug.J*X);
    BX = aug.E * X;
    gamma = X' * aug.bb;
    KX = [A, BX; sparse(s,n), G];

    aug.KX = KX;
    aug.cX = [b0; gamma];
    aug.XK = XK;
    aug.X = X;
    aug.BX = BX;
    aug.G = G;
    aug.gamma = gamma;
    aug.QX = [speye(n), sparse(n,s); sparse(n*s,n), sparse(X)];

    if build_all
        % M_K is the metric matrix for the KIOPS coordinates.
        RK = blkdiag(speye(n), R);
        MK = RK' * RK;

        % F_{M_K}(K) is computed as the Euclidean FOV of R_K*K*R_K^{-1}.
        RK_K_RKinv = RK * aug.K / RK;

        aug.QK = [speye(n), sparse(n,s); sparse(n*s,n), sparse(XK)];
        aug.RK = RK;
        aug.MK = MK;
        aug.RK_K_RKinv = RK_K_RKinv;
    end
end

end
