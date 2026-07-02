function family = make_rhs_family(n, s)
%MAKE_RHS_FAMILY Draw shared random data for right-hand-side sensitivity tests.

family.n = n;
family.s = s;

b0 = randn(n, 1);
family.b0 = b0 / norm(b0);

% Common direction q shared by all b_j before perturbations are added.
q = randn(n, 1);
family.q = q / norm(q);

family.r = zeros(n, s);
for j = 1:s
    % Independent perturbation directions r_j are later scaled by delta.
    rj = randn(n, 1);
    family.r(:,j) = rj / norm(rj);
end

end
