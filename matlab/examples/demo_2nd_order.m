% Demo: 2nd-order LCCDE via direct and unilateral methods
% Equation: y[n] - 0.5 y[n-1] - 0.25 y[n-2] = x[n]
% a = [1, -0.5, -0.25], b = [1]

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));

n = -10:60;                    % index support
x = dt_signals('step', n);     % unit step input

a = [1, -0.5, -0.25];
b = [1];

% Initial conditions: y[-1] = 0, y[-2] = 0 (zero ICs)
ic = struct('type', 'past', 'values', [0 0]);

% Direct method
sol_d = lccde_direct(a, b, n, x, ic);

% Unilateral (impulse-response + convolution)
sol_u = lccde_unilateral(a, b, n, x, ic);

figure('Name','LCCDE 2nd-Order Demo');
subplot(3,1,1);
stem(n, x, 'filled'); grid on; title('Input x[n] (unit step)'); xlabel('n'); ylabel('x[n]');

subplot(3,1,2);
stem(n, sol_d.y_h, 'r'); hold on; stem(n, sol_d.y_p, 'b'); grid on;
legend('y_h (direct)','y_p (direct)'); title('Direct method components'); xlabel('n'); ylabel('amplitude');

subplot(3,1,3);
stem(n, sol_d.y, 'k'); hold on; stem(n, sol_u.y, 'g.'); grid on;
legend('total (direct)','total (unilateral)'); title('Total outputs comparison'); xlabel('n'); ylabel('y[n]');
