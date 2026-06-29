function est = converg_bestpoly_bnd(boundary, beta, norm_ref, mvals, c_cr)
%CONVERG_BESTPOLY_BND  FOV convergence estimates using best polynomial 
% samples.

z = boundary(:);
est = zeros(size(mvals));

opts.InitialN = 1000;
opts.CheckN = 10000;
opts.AdaptiveIters = 1;
opts.AddWorst = 30;
opts.NumDirections = 32;
opts.LPRefineSteps = 2;
opts.Display = "none";
opts.Verbose = false;

for k = 1:length(mvals)
    m = mvals(k);
    out = bestpoly_bnd_fov(z, m-1, opts);
    est(k) = c_cr * beta * out.upper / norm_ref;
end

est(~isfinite(est)) = realmax;

end
