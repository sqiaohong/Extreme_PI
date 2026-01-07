function y_smooth = spline_random_smooth(x, y, n_repeat)


if nargin < 3, n_repeat = 10; end
x = x(:); y = y(:);
ny = numel(x);
y_all = zeros(ny, n_repeat);

for k = 1:n_repeat
   
    knot_spacing = 15 + (30-15)*rand;
    n_knots = max(3, round((x(end)-x(1))/knot_spacing));
    knots = linspace(x(1), x(end), n_knots+2);
    knots = knots(2:end-1);  

    pp = csaps(x, y, 1);     
    y_all(:,k) = fnval(pp, x);
end

y_smooth = mean(y_all, 2);
end
