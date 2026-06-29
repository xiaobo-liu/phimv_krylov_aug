function [V, H, m_done, beta] = arnoldi_mgs(A, b, m, tol)
%ARNOLDI_MGS  Arnoldi iteration via modified Gram-Schmidt with 
% reorthogonalization.

if nargin < 4 || isempty(tol)
    tol = 0;
end

N = size(A,1);
beta = norm(b);

V = zeros(N, m+1);
H = zeros(m+1, m);
m_done = m;

if beta == 0
    m_done = 0;
    return
end

V(:,1) = b / beta;

for j = 1:m
    w = A * V(:,j);

    % Modified Gram-Schmidt.
    for i = 1:j
        H(i,j) = V(:,i)' * w;
        w = w - H(i,j) * V(:,i);
    end

    % A second MGS pass reduces loss of orthogonality in finite precision.
    for i = 1:j
        h = V(:,i)' * w;
        H(i,j) = H(i,j) + h;
        w = w - h * V(:,i);
    end

    H(j+1,j) = norm(w);
    if H(j+1,j) <= tol
        m_done = j;
        break
    end

    V(:,j+1) = w / H(j+1,j);
end

end
