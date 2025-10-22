% classify_discrete_time_signal.m
% Interactive classifier for discrete-time signals x[n]
% - Enter x[n] as a function of n using element-wise operators (.^, .*, ./)
% - Supports helper functions: u(n) = unit step, delta(n) = Kronecker delta
% - The script estimates partial energy E_N and average power P_N numerically
%   over symmetric windows [-N, N], then classifies as ENERGY, POWER, or NEITHER
%
% Notes:
% - Numerical classification depends on a finite window. Increase N if unsure.
% - Zero signal is reported as ENERGY (convention: both energy and power are zero).
%
% Examples to try:
%   (0.5).^n .* u(n)                 % right-sided decaying exponential -> ENERGY
%   cos(0.2*pi*n)                    % periodic with finite power -> POWER
%   double(n==0)                     % delta[n] -> ENERGY (E = 1)
%   (1./sqrt(abs(n)+1)) .* (n~=0)    % diverging energy, vanishing power -> NEITHER

fprintf('Enter x[n] as a function of n (use element-wise ops .^, .*, ./).\n');
fprintf('Helpers available: u(n) for unit step, delta(n) for Kronecker delta.\n');
fprintf('Examples: (0.5).^n .* u(n), cos(0.2*pi*n), double(n==0)\n\n');

expr = input('x[n] = ', 's');

% Choose the maximum |n| to examine (window half-length)
DEFAULT_NMAX = 5000;
userN = input(sprintf('Max |n| to examine [default %d]: ', DEFAULT_NMAX));
if isempty(userN)
    Nmax = DEFAULT_NMAX;
else
    Nmax = max(10, round(abs(userN)));
end

% Build n-grid and helpers in the current workspace.
n = (-Nmax:Nmax).';             % column vector for broadcasting
idxZero = Nmax + 1;             % index where n == 0

% Helper functions (available to user expressions via eval)
u = @(k) double(k >= 0);        % unit step u[n]
delta = @(k) double(k == 0);    % Kronecker delta

% Evaluate the expression across the grid.
try
    x = eval(expr); %#ok<EVLDIR> % evaluate using current workspace (n, u, delta)
catch ME
    fprintf(2, '\nCould not evaluate expression. Error:\n  %s\n', ME.message);
    return;
end

% Ensure x is a vector of same size as n.
if isscalar(x)
    x = x + zeros(size(n));
end
if ~isequal(size(x), size(n))
    fprintf(2, '\nExpression must evaluate to a vector of size %d-by-1.\n', numel(n));
    return;
end

if any(~isfinite(x))
    warning('Expression produced Inf/NaN within the examined range; results may be unreliable.');
end

powerSamples = abs(x).^2;

% Choose a set of N values to probe convergence (spanning the window)
fractions = [0.02 0.05 0.1 0.2 0.4 0.6 0.8 1.0];
Nlist = unique(max(1, round(fractions * Nmax)));
if numel(Nlist) < 2
    Nlist = unique([max(1, round(0.5*Nmax)), Nmax]);
end

Evals = zeros(numel(Nlist), 1);
Pvals = zeros(numel(Nlist), 1);
for i = 1:numel(Nlist)
    Ni = Nlist(i);
    a = idxZero - Ni;
    b = idxZero + Ni;
    Ei = sum(powerSamples(a:b));
    Evals(i) = Ei;
    Pvals(i) = Ei / (2*Ni + 1);
end

E_last = Evals(end);
P_last = Pvals(end);
E_prev = Evals(max(1, end-1));
P_prev = Pvals(max(1, end-1));

% Numerical thresholds
SMALL = 1e-12;
ENERGY_TAIL_TOL = 5e-4;  % relative change tolerance on E_N tail
POWER_TOL = 1e-3;        % relative stability tolerance on P_N tail
P_MIN = 1e-8;            % minimum nonzero power to avoid classifying as power when P->0
E_CAP = 1e12;            % guard against numeric overflow

isZero = E_last < SMALL;
relTail = (E_last - E_prev) / max(E_last, 1);
energyConverged = (E_last < E_CAP) && (abs(relTail) < ENERGY_TAIL_TOL);
powerStable = (abs(P_last - P_prev) <= max(1e-9, POWER_TOL * max(P_last, 1)));
isPower = (~energyConverged) && powerStable && (P_last > P_MIN);

fprintf('\n');
if isZero
    fprintf('Classification: ENERGY (zero signal). E ≈ 0, P ≈ 0\n');
elseif energyConverged
    fprintf('Classification: ENERGY signal.\n');
    fprintf('Estimated energy E ≈ %.10g\n', E_last);
    fprintf('Estimated average power P_N at N=%d ≈ %.10g (→ 0 for energy signals)\n', Nlist(end), P_last);
elseif isPower
    fprintf('Classification: POWER signal.\n');
    fprintf('Estimated average power P ≈ %.10g\n', P_last);
    if isfinite(E_last)
        fprintf('Partial energy at N=%d: E_N ≈ %.10g (grows without bound)\n', Nlist(end), E_last);
    end
else
    fprintf('Classification: NEITHER energy nor power signal.\n');
    fprintf('Partial energy at N=%d: E_N ≈ %.10g\n', Nlist(end), E_last);
    fprintf('Average power estimate at N=%d: P_N ≈ %.10g\n', Nlist(end), P_last);
end

% Optional plotting for diagnosis
choice = input('Plot E_N and P_N vs N? [y/N]: ', 's');
if ~isempty(choice) && (choice == 'y' || choice == 'Y')
    figure('Name', 'Energy/Power Estimates', 'Color', 'w');
    subplot(2,1,1);
    plot(Nlist, Evals, 'o-','LineWidth',1.2); grid on;
    xlabel('N'); ylabel('E_N = \Sigma_{n=-N}^{N} |x[n]|^2');
    title('Partial energy E_N');
    subplot(2,1,2);
    plot(Nlist, Pvals, 'o-','LineWidth',1.2); grid on;
    xlabel('N'); ylabel('P_N = E_N/(2N+1)');
    title('Average power estimate P_N');
end
