function [U, family] = rhs_vectors(n, s, beta, delta, family)
%RHS_VECTORS Right-hand-side vectors U = [b0,b1,...,bs].
%
%   U = RHS_VECTORS(n,s,beta,delta) draws fresh random data and returns
%   b_j = beta*(q + delta*r_j)/norm(q + delta*r_j).
%
%   [U,family] = RHS_VECTORS(...,family) reuses the supplied random data,
%   which is useful when varying beta or delta independently.

if nargin < 5 || isempty(family)
    family = make_rhs_family(n, s);
else
    validate_rhs_family(family, n, s);
end

U = zeros(n, s+1);
U(:,1) = family.b0;

for j = 1:s
    % beta sets the common norm, while delta sets the angular perturbation.
    bj = family.q + delta*family.r(:,j);
    bj_norm = norm(bj);
    U(:,j+1) = beta * bj / bj_norm;
end

end

function validate_rhs_family(family, n, s)
%VALIDATE_RHS_FAMILY Check that reusable right-hand-side data matches this experiment.

if ~isstruct(family) || ~isfield(family, 'n') || ~isfield(family, 's') || ...
        ~isfield(family, 'b0') || ~isfield(family, 'q') || ...
        ~isfield(family, 'r')
    % Required format:
    %   family.n   scalar dimension n
    %   family.s   scalar number of perturbation directions
    %   family.b0  n-by-1 starting vector
    %   family.q   n-by-1 common direction
    %   family.r   n-by-s matrix with perturbation directions in columns
    error('The right-hand-side family must be a struct returned by make_rhs_family.');
end

if family.n ~= n || family.s ~= s
    error('The right-hand-side family has size n=%d, s=%d, but n=%d, s=%d was requested.', ...
        family.n, family.s, n, s);
end

end
