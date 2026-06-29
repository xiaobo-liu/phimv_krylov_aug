function z = fov_boundary(A, opt)
%FOV_BOUNDARY Sample the boundary of the field of values F(A).

num_theta = opt.num_theta;
theta = linspace(0, 2*pi, num_theta+1);
theta(end) = [];

A = full(A);
z = zeros(1, num_theta);

for k = 1:num_theta
    t = theta(k);

    % The largest eigenvector of this Hermitian part gives the boundary
    % point where the supporting line touches the field of values.
    H = (exp(-1i*t)*A + exp(1i*t)*A')/2;
    H = (H + H')/2;
    [V,D] = eig(H);
    [~,idx] = max(real(diag(D)));
    v = V(:,idx);

    z(k) = v' * A * v;
end

end
