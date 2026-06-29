function [A, label] = testmat(matname, n, matrix_scale)
%TESTMAT  Gallery matrices used in the FOV/Krylov experiments.

switch lower(matname)
    
    case {'laplace', 'poisson'}
        laplace_grid = 7;
        A0 = gallery('poisson', laplace_grid);
        sign = -1;
        label = sprintf('-scaled gallery(''poisson'', %d)', laplace_grid);
        
    case 'kms'
        rho = 0.88;
        A0 = gallery('kms', n, rho);
        sign = -1;
        label = sprintf('-scaled gallery(''kms'', %d, %.2f)', n, rho);
       
    case 'grcar'
        A0 = gallery('grcar', n) - 2*eye(n);
        sign = 1;
        label = sprintf('shifted/scaled gallery(''grcar'', %d)', n);

    case 'dorr'
        A0 = gallery('dorr', n);
        sign = -1;
        label = sprintf('-scaled gallery(''dorr'', %d)', n);

    case 'triw'
        A0 = gallery('triw', n) - 1.5*eye(n);
        sign = 1;
        label = sprintf('shifted/scaled gallery(''triw'', %d)', n);

    otherwise
        error('Unknown test matrix "%s".', matname);
end

% The scale keeps expm(K) well behaved while preserving FOV shape.
A = sparse(sign * matrix_scale * A0 / norm(full(A0), 2));

end
