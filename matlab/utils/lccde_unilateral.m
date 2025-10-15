function sol = lccde_unilateral(a, b, n, x, ic)
%LCCDE_UNILATERAL Solve LCCDE via unilateral z-transform approach (numeric).
%   sol = LCCDE_UNILATERAL(a, b, n, x, ic)
%   Computes zero-state response via convolution with impulse response h[n],
%   and zero-input response via homogeneous recursion with initial conditions.
%
%   Inputs/Outputs: same conventions as lccde_direct.m
%
arguments
    a (1,:) double
    b (1,:) double
    n (1,:) double
    x (1,:) double
    ic struct
end

N = numel(a) - 1; %#ok<NASGU>
assert(a(1) ~= 0, 'a(1) (a0) must be nonzero.');
assert(isvector(n) && isvector(x) && numel(n)==numel(x), 'n and x must match.');

% 1) Impulse response h via direct solver with delta input, zero ICs
x_delta = double(n==0);
ics_zero = struct('type', 'past', 'values', zeros(1, max(numel(a)-1,0)));
h_sol = lccde_direct(a, b, n, x_delta, ics_zero);
h = h_sol.y_p; % zero-state response to delta is impulse response

% 2) Zero-state response to arbitrary input via convolution
[y_conv, ny] = dt_conv(x, n, h, n);
% Sample onto requested n support
% Build lookup from ny to indices
ny_min = ny(1); ny_max = ny(end);
ys_on_n = zeros(size(n));
for i = 1:numel(n)
    if n(i) < ny_min || n(i) > ny_max
        ys_on_n(i) = 0;
    else
        ys_on_n(i) = y_conv(n(i) - ny_min + 1);
    end
end

% 3) Zero-input response via homogeneous recursion (direct solver with x=0)
x_zero = zeros(size(x));
hom_sol = lccde_direct(a, b, n, x_zero, ic);
y_h = hom_sol.y_h; % zero-input part

% 4) Total
y_total = ys_on_n + y_h;

sol = struct('y', y_total, 'y_h', y_h, 'y_p', ys_on_n, 'h', h, 'n', n, 'ny_conv', ny);
end
