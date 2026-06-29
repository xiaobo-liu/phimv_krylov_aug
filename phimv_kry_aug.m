function [y, est] = phimv_kry_aug(A, U, s, m, tol, basis)
%PHIMV_KRY_AUG Compute phi-combination action using a chosen basis.
%
%   [y,est] = phimv_kry_aug(A,U,s,m,tol,basis) computes
%               y = sum_{i=0}^{s} phi_i(A) b_i,
%   where U is an n-by-(s+1) matrix with columns b_0,...,b_s.
%
%   The basis argument chooses the augmented realization used by Arnoldi:
%       'kiops' or 'k'   Arnoldi on the KIOPS coordinate matrix K.
%       'block' or 'w'   Arnoldi on the large block matrix W.
%       'orth'  or 'kx'  Arnoldi on the orthogonally compressed matrix K_X.
%
%   If basis is omitted, the default is 'kiops'.  For convenience,
%   phimv_kry_aug(A,U,s,m,basis) is also accepted.

if nargin < 5
    tol = [];
end
if nargin < 6
    basis = 'kiops';
end

% Allow phimv_kry_aug(A,U,s,m,'orth') with the default tolerance.
if ischar(tol) || isa(tol, 'string')
    basis = char(tol);
    tol = [];
end

if isempty(tol)
    tol = eps('double')/2;
end
if isempty(basis)
    basis = 'kiops';
end
basis = lower(char(basis));

[n, n2] = size(A);
if n ~= n2
    error('A must be square.');
end

[nU, cols] = size(U);
if s < 1
    error('s must be no smaller than 1.');
end
if nU ~= n || s ~= cols - 1
    error('U must be n-by-(s+1) with same n as A.');
end

aug = build_aug_mat(A, U, basis);

switch basis
    case {'kiops', 'k'}
        M = aug.K;
        b = aug.c;
    case {'block', 'w'}
        M = aug.W;
        b = aug.b;
    case {'orth', 'kx'}
        M = aug.KX;
        b = aug.cX;
    otherwise
        error('Unknown basis "%s". Use "kiops"/"k", "block"/"w", or "orth"/"kx".', basis);
end

beta = norm(b,2);
if beta == 0
    y = zeros(n,1);
    est = 0;
    return
end

% Modified Gram-Schmidt Arnoldi.
[V, H, m, beta] = arnoldi_mgs(M, b, m, tol);

% exp(H_m)e_1 gives the coefficient vector in the Arnoldi basis.
e1 = zeros(m,1);
e1(1) = 1;
g = expm(H(1:m,1:m)) * e1;

% Only the first n entries represent the desired linear combination.
V_top = V(1:n, 1:m);
y = beta * (V_top * g);

% Error estimate beta * h_{m+1,m} * |g_m|
est = beta * H(m+1,m) * abs(g(m));
end
