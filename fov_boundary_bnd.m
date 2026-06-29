function z = fov_boundary_bnd(points, radius, num_theta)
%FOV_BOUNDARY_BND Boundary of the FOV bound conv(points) + disk(0,radius).

theta = linspace(0, 2*pi, num_theta+1);
theta(end) = [];

points = points(:);
h = zeros(size(theta));

for k = 1:num_theta
    % Support functions add under Minkowski sums.
    h(k) = max(real(exp(-1i*theta(k))*points)) + radius;
end

dtheta = theta(2) - theta(1);
hp = (circshift(h,-1) - circshift(h,1)) / (2*dtheta);

% Smooth support-function parametrization of the boundary.
z = exp(1i*theta) .* (h + 1i*hp);

end
