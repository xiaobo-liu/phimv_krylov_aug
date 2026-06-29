function out = bestpoly_bnd_fov(z, deg, opts)
%BESTPOLY_BND_FOV Adaptive sampled best-polynomial bound on a FOV boundary.
%
%   out = bestpoly_bnd_fov(z,deg,opts) approximates the degree-deg
%   minimax error for exp on the polygonal boundary z.  The experiments use
%   out.upper, the maximum error on a dense validation sample.

    if nargin < 3
        opts = struct();
    end

    % Sampling, angular discretization, and linear-programming defaults.
    opts = setDefault(opts, "InitialN", 1200);
    opts = setDefault(opts, "CheckN", 12000);
    opts = setDefault(opts, "AdaptiveIters", 5);
    opts = setDefault(opts, "AddWorst", 80);
    opts = setDefault(opts, "NumDirections", 64);
    opts = setDefault(opts, "LPRefineSteps", 4);
    opts = setDefault(opts, "TolUnique", 1e-13);
    opts = setDefault(opts, "Display", "none");
    opts = setDefault(opts, "Algorithm", "dual-simplex");
    opts = setDefault(opts, "ConstraintTolerance", 1e-10);
    opts = setDefault(opts, "OptimalityTolerance", 1e-10);
    opts = setDefault(opts, "Verbose", true);

    z = closePolygon(z(:));

    % Center and scale the polygon before building the Arnoldi basis.
    center = 0.5*(min(real(z)) + max(real(z))) ...
           + 0.5i*(min(imag(z)) + max(imag(z)));
    rho = max(abs(z - center));
    if rho == 0
        rho = 1;
    end

    zgrid = samplePolygonTotal(z, opts.InitialN);
    % Always include the vertices, since FOV corners can control the error.
    zgrid = uniqueComplexTol([zgrid; z(1:end-1)], opts.TolUnique);

    % Direction angles approximate the complex modulus by supporting lines.
    theta = 2*pi*(0:opts.NumDirections-1).' / opts.NumDirections;

    hist = struct([]);

    for it = 1:opts.AdaptiveIters
        % Solve on the current grid, then validate on a denser sample.
        [coef, basis, fscale, gridErr, lpTau, exitflag, output] = ...
            solveGridLPRefined(zgrid, deg, center, rho, theta, opts);

        pfun = @(w) evalArnoldiPoly(w, coef, basis, center, rho, fscale);

        zcheck = samplePolygonTotal(z, opts.CheckN);
        err = abs(exp(zcheck) - pfun(zcheck));
        [upper, imax] = max(err);

        hist(it).numGrid = numel(zgrid);
        hist(it).gridErr = gridErr;
        hist(it).lpTau = lpTau;
        hist(it).upper = upper;
        hist(it).ratio = upper / max(gridErr, realmin);
        hist(it).zmax = zcheck(imax);
        hist(it).exitflag = exitflag;
        hist(it).output = output;

        if opts.Verbose
            fprintf("it %2d: grid=%5d  gridErr=%.6e  checkErr=%.6e  ratio=%.4f\n", ...
                it, numel(zgrid), gridErr, upper, upper/max(gridErr,realmin));
        end

        % Add the largest local error peaks to the next grid.
        idx = selectWorstLocalMaxima(err, opts.AddWorst);
        zgrid = uniqueComplexTol([zgrid; zcheck(idx); z(1:end-1)], opts.TolUnique);
    end

    % Re-solve after the final adaptive enrichment.
    [coef, basis, fscale, gridErr, lpTau, exitflag, output] = ...
        solveGridLPRefined(zgrid, deg, center, rho, theta, opts);

    pfun = @(w) evalArnoldiPoly(w, coef, basis, center, rho, fscale);

    zcheck = samplePolygonTotal(z, opts.CheckN);
    err = abs(exp(zcheck) - pfun(zcheck));
    [upper, imax] = max(err);

    out.coef = coef;
    out.basis = basis;
    out.fscale = fscale;
    out.center = center;
    out.rho = rho;
    out.grid = zgrid;
    out.gridErr = gridErr;
    out.lpTau = lpTau;
    out.upper = upper;
    out.zmax = zcheck(imax);
    out.history = hist;
    out.exitflag = exitflag;
    out.output = output;
    out.angularSlack = 1/cos(pi/opts.NumDirections);
    % Function handle for checking or plotting the computed polynomial.
    out.eval = pfun;

    if opts.Verbose
        fprintf("final grid error      = %.16e\n", gridErr);
        fprintf("final validation error= %.16e\n", upper);
        fprintf("angular slack factor  = %.16e\n", out.angularSlack);
    end
end

function [coef, basis, fscale, gridErr, lpTau, exitflag, output] = ...
    solveGridLPRefined(zgrid, deg, center, rho, theta, opts)
%SOLVEGRIDLPREFINED Solve the sampled problem and LP corrections.

    xi = (zgrid - center) / rho;
    % The discrete Arnoldi basis is more stable than monomials on the grid.
    [basis, P] = discreteArnoldiBasis(xi, deg);

    fraw = exp(zgrid);
    fscale = max(1, max(abs(fraw)));
    f = fraw / fscale;

    coef = P \ f;

    % Start from least squares, then reduce the maximum sampled residual.
    bestCoef = coef;
    bestErr = max(abs(P*coef - f));

    lpTau = NaN;
    exitflag = NaN;
    output = [];

    for ref = 1:opts.LPRefineSteps
        rcur = f - P*coef;
        rscale = max(abs(rcur));

        if rscale <= 10*eps(max(1,norm(f,inf)))
            break
        end

        g = rcur / rscale;

        % The LP correction is solved for the normalized residual.
        [dcoef, tau, exitflag, output] = solveCorrectionLP(P, g, theta, opts);

        if isempty(dcoef) || any(~isfinite(dcoef))
            break
        end
        
        candCoef = coef + rscale*dcoef;
        candErr = max(abs(P*candCoef - f));
        
        if ~isfinite(candErr)
            break
        end

        if candErr < bestErr
            coef = candCoef;
            bestCoef = candCoef;
            bestErr = candErr;
            lpTau = tau * rscale * fscale;
        else
            break
        end
    end

    coef = bestCoef;
    gridErr = bestErr * fscale;
end

function [dcoef, tau, exitflag, output] = solveCorrectionLP(P, g, theta, opts)
%SOLVECORRECTIONLP Minimize a sampled infinity norm by linear programming.
    [N, m] = size(P);
    L = numel(theta);

    % Variables are real(dcoef), imag(dcoef), and the scalar tau.
    nvar = 2*m + 1;
    ncon = N*L;

    Aineq = zeros(ncon, nvar);
    bineq = zeros(ncon, 1);

    Pr = real(P);
    Pi = imag(P);
    gr = real(g);
    gi = imag(g);

    row0 = 0;

    for ell = 1:L
        % real(exp(-1i*theta)*(P*dcoef - g)) <= tau.
        cth = cos(theta(ell));
        sth = sin(theta(ell));

        rows = row0 + (1:N);

        Aineq(rows, 1:m)       =  cth*Pr + sth*Pi;
        Aineq(rows, m+1:2*m)   = -cth*Pi + sth*Pr;
        Aineq(rows, end)       = -1;

        bineq(rows) = cth*gr + sth*gi;

        row0 = row0 + N;
    end

    fobj = zeros(nvar,1);
    fobj(end) = 1;

    lb = -inf(nvar,1);
    ub =  inf(nvar,1);
    lb(end) = 0;

    options = optimoptions("linprog", ...
        "Algorithm", opts.Algorithm, ...
        "Display", opts.Display);

    options = trySetOption(options, "ConstraintTolerance", opts.ConstraintTolerance);
    options = trySetOption(options, "OptimalityTolerance", opts.OptimalityTolerance);

    try
        [x, tau, exitflag, output] = linprog( ...
            fobj, Aineq, bineq, [], [], lb, ub, options);
    catch ME
        x = [];
        tau = NaN;
        exitflag = -1000;

        output = struct();
        output.message = "linprog failed with error: " + string(ME.message);
    end

    badSolution = isempty(x) || numel(x) < nvar || exitflag <= 0;

    if badSolution
        % Try a second standard algorithm before falling back to least squares.
        try
            options2 = optimoptions("linprog", ...
                "Algorithm", "interior-point", ...
                "Display", opts.Display);

            options2 = trySetOption(options2, "ConstraintTolerance", opts.ConstraintTolerance);
            options2 = trySetOption(options2, "OptimalityTolerance", opts.OptimalityTolerance);

            [x2, tau2, exitflag2, output2] = linprog( ...
                fobj, Aineq, bineq, [], [], lb, ub, options2);

            if ~isempty(x2) && numel(x2) >= nvar && exitflag2 > 0
                x = x2;
                tau = tau2;
                exitflag = exitflag2;
                output = output2;
                badSolution = false;
            end
        catch

        end
    end


    if badSolution
        % Keep the outer refinement usable even when linprog fails.
        dcoef = P \ g;
        residual = P*dcoef - g;
        tau = max(abs(residual));
        exitflag = -999;

        output = struct();
        output.message = "linprog failed; used least-squares correction fallback.";

        return
    end

    dcoef = x(1:m) + 1i*x(m+1:2*m);
end


function [basis, Q] = discreteArnoldiBasis(xi, deg)
%DISCRETEARNOLDIBASIS Orthonormal polynomial basis on sampled points.

    xi = xi(:);
    N = numel(xi);

    if N <= deg+1
        error("Need more sample points than polynomial coefficients.");
    end

    Q = zeros(N, deg+1);
    H = zeros(deg+1, deg);

    beta = sqrt(N);
    Q(:,1) = ones(N,1) / beta;

    for k = 1:deg
        v = xi .* Q(:,k);

        % Two MGS passes reduce loss of orthogonality on clustered samples.
        for pass = 1:2
            h = Q(:,1:k)' * v;
            H(1:k,k) = H(1:k,k) + h;
            v = v - Q(:,1:k)*h;
        end

        H(k+1,k) = norm(v);

        if H(k+1,k) < 1e-14
            error("Arnoldi breakdown: increase grid size or reduce degree.");
        end

        Q(:,k+1) = v / H(k+1,k);
    end

    basis.H = H;
    basis.beta = beta;
    basis.deg = deg;
end

function y = evalArnoldiPoly(z, coef, basis, center, rho, fscale)
%EVALARNOLDIPOLY Evaluate the scaled polynomial in physical coordinates.

    sz = size(z);
    xi = (z(:)-center)/rho;

    Q = evalArnoldiBasis(xi, basis);
    y = fscale * (Q*coef);

    y = reshape(y, sz);
end

function Q = evalArnoldiBasis(xi, basis)
%EVALARNOLDIBASIS Evaluate the Arnoldi polynomial basis.

    xi = xi(:);
    deg = basis.deg;
    H = basis.H;

    Q = zeros(numel(xi), deg+1);
    Q(:,1) = 1 / basis.beta;

    for k = 1:deg
        v = xi .* Q(:,k);

        for j = 1:k
            v = v - H(j,k)*Q(:,j);
        end

        Q(:,k+1) = v / H(k+1,k);
    end
end

function zsam = samplePolygonTotal(z, Ntotal)
%SAMPLEPOLYGONTOTAL Sample a closed polygon uniformly in arclength.

    z = closePolygon(z(:));

    edge = z(2:end) - z(1:end-1);
    len = abs(edge);
    total = sum(len);

    if total == 0
        zsam = z(1);
        return
    end

    cumlen = [0; cumsum(len)];

    t = (0:Ntotal-1).' / Ntotal * total;

    idx = discretize(t, cumlen);
    idx(isnan(idx)) = numel(len);
    idx = min(idx, numel(len));

    local = (t - cumlen(idx)) ./ len(idx);
    zsam = z(idx) + local .* edge(idx);
end

function idx = selectWorstLocalMaxima(err, k)
%SELECTWORSTLOCALMAXIMA Select largest local error peaks for refinement.

    err = err(:);
    n = numel(err);

    idx = find(err >= err([n,1:n-1]) & err >= err([2:n,1]));

    if isempty(idx)
        [~, idx] = sort(err, "descend");
    else
        [~, ord] = sort(err(idx), "descend");
        idx = idx(ord);
    end

    idx = idx(1:min(k,numel(idx)));
end

function z = closePolygon(z)
%CLOSEPOLYGON Append the first point if the boundary is not closed.

    if isempty(z)
        error("Boundary vector is empty.");
    end

    scale = max(1, norm(z,inf));

    if abs(z(end)-z(1)) > 100*eps(scale)
        z = [z; z(1)];
    end
end

function zuniq = uniqueComplexTol(z, tol)
%UNIQUECOMPLEXTOL Remove nearly duplicate complex points.

    z = z(:);
    scale = max(1, max(abs(z)));
    h = tol * scale;

    key = [round(real(z)/h), round(imag(z)/h)];
    [~, ia] = unique(key, "rows", "stable");

    zuniq = z(ia);
end

function opts = setDefault(opts, name, value)
%SETDEFAULT Set a structure field when it is absent.

    if ~isfield(opts, name)
        opts.(name) = value;
    end
end

function options = trySetOption(options, name, value)
%TRYSETOPTION Set an optimizer option when this MATLAB release supports it.

    try
        options = optimoptions(options, name, value);
    catch

    end
end
