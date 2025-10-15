function x = dt_signals(kind, n, varargin)
%DT_SIGNALS Generate common discrete-time signals.
%   x = DT_SIGNALS(kind, n, params...) supports:
%     kind:
%       'delta'            : unit impulse δ[n]
%       'step'             : unit step u[n]
%       'ramp'             : discrete ramp n*u[n]
%       'exp'              : exponential a.^n * u[n] (param: a)
%       'sin'              : sin(Ω0*n + φ) * u[n] (params: Omega0, phi)
%       'cos'              : cos(Ω0*n + φ) * u[n] (params: Omega0, phi)
%       'custom-pw'        : piecewise from segments {n0,n1, valFcnHandle}
%
%   n is integer vector (e.g., -50:50).
%
%   Examples:
%     n = -10:40; x = dt_signals('exp', n, 0.8);
%     x = dt_signals('sin', n, pi/4, 0);
%     x = dt_signals('delta', n);
%
%   This utility applies u[n]=1 for n>=0 when appropriate.

arguments
    kind (1,:) char
    n (1,:) double
end

switch lower(kind)
    case 'delta'
        x = double(n == 0);

    case 'step'
        x = double(n >= 0);

    case 'ramp'
        x = n .* double(n >= 0);

    case 'exp'
        a = varargin{1};
        x = (a .^ n) .* double(n >= 0);

    case {'sin','sine'}
        Omega0 = varargin{1};
        phi    = 0; if numel(varargin) >= 2, phi = varargin{2}; end
        x = sin(Omega0 .* n + phi) .* double(n >= 0);

    case 'cos'
        Omega0 = varargin{1};
        phi    = 0; if numel(varargin) >= 2, phi = varargin{2}; end
        x = cos(Omega0 .* n + phi) .* double(n >= 0);

    case 'custom-pw'
        % varargin: cell array of segments {n0, n1, f(n)}
        assert(~isempty(varargin), 'Provide segments for custom-pw');
        segments = varargin{1};
        x = zeros(size(n));
        for k = 1:numel(segments)
            seg = segments{k};
            n0 = seg{1}; n1 = seg{2}; f = seg{3};
            idx = (n >= n0) & (n <= n1);
            x(idx) = f(n(idx));
        end

    otherwise
        error('Unsupported kind: %s', kind);
end
