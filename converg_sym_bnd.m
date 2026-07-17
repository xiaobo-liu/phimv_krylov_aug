function est = converg_sym_bnd(boundary, start_norm, norm_ref, mvals, c_cr, s)
%CONVERG_SYM_BND  Explicit Chebyshev bound for symmetric A (Proposition
%                 5.5. in the paper)

lmax = max(abs(boundary));
est = zeros(size(mvals));

for k = 1:length(mvals)
    m = mvals(k);
    out = cheb_bnd_fov_sym(lmax, s, m);
    est(k) = c_cr * start_norm* out.bound / norm_ref;
end

est(~isfinite(est)) = realmax;

end
