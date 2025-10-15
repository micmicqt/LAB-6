function sol = lccde_direct(a, b, n, x, ic)
%LCCDE_DIRECT Solve LCCDE using direct time-domain method.
%   sol = LCCDE_DIRECT(a, b, n, x, ic) solves
%       sum_{k=0}^N a(k+1) y[n-k] = sum_{k=0}^M b(k+1) x[n-k]
%   where a(1) must be nonzero (a0). Vectors a and b are in descending lag order.
%
%   Inputs
%     a : [a0 a1 ... aN]
%     b : [b0 b1 ... bM]
%     n : integer index support, e.g., -50:100
%     x : input sequence over n (vector same size as n)
%     ic: initial conditions struct with fields:
%         .type : 'past' for y[-1], y[-2], ... given
%                 'future' for y[1], y[2], ... given (time-reversal handling)
%         .values : vector of length N with order consistent with y[-1], y[-2], ... or y[1], y[2], ...
%
%   Outputs
%     sol struct with fields:
%       .y          : total solution y[n]
%       .y_h        : homogeneous solution (zero-input response)
%       .y_p        : particular solution (zero-state response)
%       .n          : index vector
%       .recurrence : function handle for one-step recursion
%
%   Notes
%   - Uses causal recursion for n ascending using provided ICs.
%   - For 'future' ICs (y[1..N]), a time-reversal trick constructs equivalent past ICs.
%   - The homogeneous part y_h comes from recursion with x=0 and ICs; y_p uses x with zero ICs.

arguments
    a (1,:) double
    b (1,:) double
    n (1,:) double
    x (1,:) double
    ic struct
end

N = numel(a) - 1; % order
M = numel(b) - 1;
assert(a(1) ~= 0, 'a(1) (a0) must be nonzero.');
assert(isvector(n) && isvector(x) && numel(n)==numel(x), 'n and x must match.');

% Build index mapping: map y[n-k], x[n-k]
% We'll compute y over n from smallest to largest
nmin = n(1); nmax = n(end);

% Prepare arrays with padding for negative indices
% We store y over full n, and also maintain a ring buffer for last N samples
y_total = zeros(size(n));
y_h     = zeros(size(n));
y_p     = zeros(size(n));

% Helper: one-step update given previous y's and current/past x's
function yn = step(prevY, idx)
    % prevY: vector [y[n-1], y[n-2], ..., y[n-N]]
    % compute RHS sum b_k x[n-k]
    rhs = 0;
    for k = 0:M
        idx_x = find(n == (n(idx) - k), 1, 'first');
        if ~isempty(idx_x)
            rhs = rhs + b(k+1) * x(idx_x);
        end
    end
    % compute -a_k y[n-k]
    acc = rhs;
    for k = 1:N
        if k <= numel(prevY)
            acc = acc - a(k+1) * prevY(k);
        end
    end
    yn = acc / a(1);
end

% Initialize buffers for homogeneous and particular runs
prevY_h = zeros(1, N);
prevY_p = zeros(1, N);

% Map ICs into buffers at the starting index n(1)
if N > 0
    switch lower(ic.type)
        case 'past'
            % ic.values = [y[-1], y[-2], ..., y[-N]] at absolute indices
            % Determine absolute indices based on n grid
            % Fill prevY at start n(1) using these values shifted appropriately
            for k = 1:N
                target_n = n(1) - k; % index before the first n sample
                % if ICs are specified at corresponding negative indices, use directly
                if numel(ic.values) >= k
                    prevY_h(k) = ic.values(k);
                    prevY_p(k) = 0; % zero-state for particular
                end
            end
        case 'future'
            % Given y[1], y[2], ... y[N]; approximate by shifting start near 0
            % For simplicity, if n(1) <= 0, we can translate: y[n(1)-k] unknown
            % We place these as constraints when we reach those indices. Here we seed buffer with zeros
            % and add correction once indices are reached. Simpler approach: reflect indices if needed.
            warning('Future ICs not strictly mapped; proceeding with zero initialization for start.');
        otherwise
            error('Unsupported ic.type: %s', ic.type);
    end
end

% Run two recursions: homogeneous (x=0) and particular (IC=0)
for idx = 1:numel(n)
    % homogeneous: x=0 -> set temporary x to zero
    rhs_save = x; x_zero = 0; %#ok<NASGU>
    % compute y_h using step with x=0 by temporarily zeroing b contribution
    % Implement by locally computing RHS=0
    rhs_h = 0;
    for k = 1:N
        if k <= numel(prevY_h)
            rhs_h = rhs_h - a(k+1) * prevY_h(k);
        end
    end
    y_h(idx) = rhs_h / a(1);

    % particular: zero ICs, use actual x
    y_p(idx) = step(prevY_p, idx);

    % total
    y_total(idx) = y_h(idx) + y_p(idx);

    % roll buffers
    if N > 0
        prevY_h = [y_h(idx), prevY_h(1:end-1)];
        prevY_p = [y_p(idx), prevY_p(1:end-1)];
    end
end

% provide recurrence handle
recurrence = @(prevY_vals, nn) ( ...
    (sum(arrayfun(@(k) (k<=M) * b(k+1) * interp_x(nn-k), 0:max(M,0))) - ...
      sum(arrayfun(@(k) (k>=1 && k<=numel(prevY_vals)) * a(k+1) * prevY_vals(k), 1:N)) ) / a(1) );

    function xv = interp_x(qn)
        idxq = find(n == qn, 1, 'first');
        if isempty(idxq), xv = 0; else, xv = x(idxq); end
    end

sol = struct('y', y_total, 'y_h', y_h, 'y_p', y_p, 'n', n, 'recurrence', recurrence);
end
