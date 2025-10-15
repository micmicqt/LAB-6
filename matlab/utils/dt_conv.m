function [y, ny] = dt_conv(x, nx, h, nh)
%DT_CONV Discrete-time linear convolution with index tracking.
%   [y, ny] = DT_CONV(x, nx, h, nh)
%   nx, nh must be integer, unit-spaced vectors.

arguments
    x (1,:) double
    nx (1,:) double
    h (1,:) double
    nh (1,:) double
end

assert(all(diff(nx)==1), 'nx must be unit-spaced contiguous integers.');
assert(all(diff(nh)==1), 'nh must be unit-spaced contiguous integers.');

% Linear convolution
y_full = conv(x, h);
ny_start = nx(1) + nh(1);
ny_end   = nx(end) + nh(end);
ny = ny_start:ny_end;

% Sanity
assert(numel(y_full) == numel(ny), 'Index vector length mismatch.');

y = y_full; % already aligned to ny
end
