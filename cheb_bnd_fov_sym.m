function out = cheb_bnd_fov_sym(alpha, s, m)
%CHEB_BND_FOV_SYM Bound from Proposition 5.5.
%
%   out = expFOVJordanZeroBound(alpha, s, m)
%
% Computes the bound
%
%   E_{m-1}(exp, F(W)) <= sum_{k=m}^inf exp(-h) B_k(h)
%                              * (chi_s^k + chi_s^(-k)),
%
% where
%
%   h     = alpha/2,
%   r_s   = cos(pi/(s+1)),
%   R_s   = r_s + 1/2,
%   chi_s = 1 + R_s/h + sqrt((R_s/h)*(2 + R_s/h)).
%
% INPUT:
%   alpha : positive number with F(A) = [-alpha,0]
%   s     : size of the nilpotent Jordan block
%   m     : Krylov dimension, i.e., the polynomial degree is m-1
%
%
% OUTPUT:
%   out.bound         : value of the bound in Proposition 5.5
%   out.converged     : whether the tail summation criterion was met
%
% The computation uses besseli(k,h,1), i.e., the exponentially scaled
% modified Bessel function exp(-h)*B_k(h), which is exactly the factor
% appearing in the proposition for real h > 0.

    tailTol = 1e-15;
    maxK = max(m + 1000, 10000);

    h  = alpha/2;
    rs = cos(pi/(s+1));
    Rs = rs + 1/2;
    chi_s = 1 + Rs/h + sqrt((Rs/h)*(2 + Rs/h));

    [bound, converged] = besselTailBound(h, chi_s, m, tailTol, maxK);

    out = struct();
    out.bound = bound;
    out.converged = converged;
end

function [tail, converged] = besselTailBound(h, chi, m, tol, maxK)
% Sum sum_{k=m}^inf exp(-h)*I_k(h)*(chi^k + chi^(-k)).

    tail = 0;
    converged = false;
    smallCount = 0;
    logChi = log(chi);

    for k = m:maxK
        scaledIk = besseli(k, h, 1); % exp(-h)*I_k(h)
        factor   = exp(k*logChi) + exp(-k*logChi);
        term     = scaledIk * factor;

        if ~isfinite(term)
            tail = Inf;
            return
        end

        tail = tail + term;

        if term <= tol*max(1, tail)
            smallCount = smallCount + 1;
        else
            smallCount = 0;
        end

        if smallCount >= 10
            converged = true;
            return
        end
    end

    warning('expFOVJordanZeroBound:tail', 'Tail summation did not meet the requested tolerance before maxK.');
end
