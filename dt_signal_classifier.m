% Discrete-Time Signal Classifier (Energy / Power / Neither)
% 
% Usage:
%   Run this script, then enter x[n] as a function of n when prompted.
%   The expression should be in terms of the index variable `n` and may use
%   the helper functions:
%     - u(n): unit step (1 for n >= 0, else 0)
%     - d(n): Kronecker delta (1 for n == 0, else 0)
% 
%   Examples to try:
%     (-0.5).^n .* u(n)
%     sin(0.2*pi*n)
%     (0.9).^abs(n)
%     d(n-3)
%     (cos(0.1*pi*n)).*(u(n) - u(n-50))
% 
% Notes:
%   - The script auto-vectorizes ^, *, / to .^, .*, ./ if needed.
%   - Classification is based on numerical tests over growing windows.
%   - Energy: sum |x[n]|^2 converges (average power tends to 0 as window grows).
%   - Power: average of |x[n]|^2 over symmetric windows converges to > 0.
%   - Neither: does not match the above within numerical tolerance.

%% Prompt for the expression (menu + optional custom)
fprintf('Discrete-Time Signal Classifier (Energy / Power / Neither)\n');
expressionString = promptExpressionViaMenu();

%% Window sizes to probe (symmetric: -N:N)
windowSizes = [64, 128, 256, 512, 1024, 2048];

%% Measure energy and average power across windows
[energyValues, averagePowerValues] = computeEnergyAndPowerMeasures(expressionString, windowSizes);

%% Classify
classificationResult = classifyFromMeasures(energyValues, averagePowerValues, windowSizes);

%% Report
fprintf('\nClassification: %s\n', classificationResult);
if strcmp(classificationResult, 'Energy signal')
    estimatedEnergy = energyValues(end);
    fprintf('Estimated energy (approx.): %.6g\n', estimatedEnergy);
elseif strcmp(classificationResult, 'Power signal')
    tailCount = min(3, numel(averagePowerValues));
    estimatedPower = mean(averagePowerValues(end-tailCount+1:end));
    fprintf('Estimated power (approx.): %.6g\n', estimatedPower);
else
    tailCount = min(3, numel(averagePowerValues));
    tailPower = mean(averagePowerValues(end-tailCount+1:end));
    fprintf('Average power over largest window (approx.): %.6g\n', tailPower);
end

%% --- Local functions ---
function [energyValues, averagePowerValues] = computeEnergyAndPowerMeasures(expressionString, windowSizes)
    numWindows = numel(windowSizes);
    energyValues = zeros(1, numWindows);
    averagePowerValues = zeros(1, numWindows);
    
    for idx = 1:numWindows
        halfWidth = windowSizes(idx);
        discreteTimeIndexVector = -halfWidth:halfWidth; % symmetric window
        signalValues = evaluateUserExpression(expressionString, discreteTimeIndexVector);
        squaredMagnitude = abs(signalValues).^2;
        squaredMagnitude(~isfinite(squaredMagnitude)) = 0; % ignore NaN/Inf
        totalEnergyOverWindow = sum(squaredMagnitude);
        energyValues(idx) = totalEnergyOverWindow;
        averagePowerValues(idx) = totalEnergyOverWindow / numel(discreteTimeIndexVector);
    end
end

function signalValues = evaluateUserExpression(expressionString, nVector)
    expressionVectorized = vectorizeExpression(expressionString);
    n = nVector; %#ok<NASGU>
    try
        signalValues = eval(expressionVectorized);
    catch evalError
        error('Failed to evaluate the expression. Vectorized form: %s\nError: %s', ...
              expressionVectorized, evalError.message);
    end
    if isscalar(signalValues)
        signalValues = signalValues .* ones(size(nVector));
    end
    if ~isvector(signalValues)
        error('The evaluated expression did not produce a 1-D vector.');
    end
    signalValues = signalValues(:).'; % force row vector
    if numel(signalValues) ~= numel(nVector)
        error('The expression produced a vector of length %d, expected %d to match n.', ...
              numel(signalValues), numel(nVector));
    end
end

function expressionVectorized = vectorizeExpression(expressionString)
    % Insert element-wise operators when the dot is missing.
    expressionVectorized = regexprep(expressionString, '(?<!\.)\^', '.^');
    expressionVectorized = regexprep(expressionVectorized, '(?<!\.)\*', '.*');
    expressionVectorized = regexprep(expressionVectorized, '(?<!\.)/', './');
end

function classificationResult = classifyFromMeasures(energyValues, averagePowerValues, windowSizes)
    smallNumber = 1e-12;
    numWindows = numel(windowSizes);
    
    if all(averagePowerValues < 100*smallNumber)
        classificationResult = 'Energy signal';
        return;
    end
    
    if numWindows < 2
        classificationResult = 'Undetermined (insufficient windows)';
        return;
    end
    
    lastN = windowSizes(end);
    prevN = windowSizes(end-1);
    lastP = averagePowerValues(end);
    prevP = averagePowerValues(end-1);
    
    predictedEnergyLastP = prevP * ((2*prevN + 1) / (2*lastN + 1));
    energyConsistencyError = abs(lastP - predictedEnergyLastP) / max(abs(lastP), smallNumber);
    powerRelativeChange = abs(lastP - prevP) / max(abs(lastP), smallNumber);
    
    if numWindows >= 3
        tailPowers = averagePowerValues(end-2:end);
        tailPowerSpread = max(tailPowers) - min(tailPowers);
        tailPowerMean = mean(tailPowers);
        powerStability = tailPowerSpread / max(abs(tailPowerMean), smallNumber);
    else
        tailPowerMean = (lastP + prevP)/2;
        powerStability = powerRelativeChange;
    end
    
    energyConsistencyThreshold = 0.05;  % close to 1/(2N+1) scaling across last step
    powerStabilityThreshold  = 0.01;    % nearly constant across last windows
    
    if (lastP < prevP) && (energyConsistencyError <= energyConsistencyThreshold)
        classificationResult = 'Energy signal';
        return;
    end
    
    if (tailPowerMean > 100*smallNumber) && (powerStability <= powerStabilityThreshold)
        classificationResult = 'Power signal';
        return;
    end
    
    classificationResult = 'Neither energy nor power signal';
end

function y = u(n)
    y = double(n >= 0);
end

function y = d(n)
    y = double(n == 0);
end

function expressionString = promptExpressionViaMenu()
    % Presents a numbered menu with common signals (no copy/paste needed).
    % Returns a MATLAB expression string in terms of n.
    while true
        fprintf('\nSelect a signal x[n] by number (or choose custom):\n');
        fprintf('  1) (-a).^n .* u(n)              (right-sided decaying/growing exponential)\n');
        fprintf('  2) sin(omega0 * n)              (discrete-time sinusoid)\n');
        fprintf('  3) (r).^abs(n)                  (two-sided exponential)\n');
        fprintf('  4) d(n - k)                     (Kronecker delta shifted by k)\n');
        fprintf('  5) cos(omega0*n).*(u(n)-u(n-N)) (finite-length cosine burst)\n');
        fprintf('  6) C                            (constant sequence)\n');
        fprintf('  7) A.*(r).^n .* u(n)            (scaled right-sided exponential)\n');
        fprintf('  8) Enter a custom expression    (type manually)\n');
        choice = input('Select option [1-8]: ');
        if isempty(choice) || ~isscalar(choice) || ~isfinite(choice)
            fprintf('Invalid selection. Try again.\n');
            continue;
        end
        choice = round(choice);
        if choice < 1 || choice > 8
            fprintf('Selection out of range. Try again.\n');
            continue;
        end
        break;
    end

    switch choice
        case 1
            a = input('Enter a (real number, e.g., 0.5): ');
            validateNumeric(a);
            expressionString = sprintf('(%s).^n .* u(n)', numToStr(a));
        case 2
            omega0 = input('Enter omega0 in radians (e.g., 0.2*pi): ');
            validateNumeric(omega0);
            expressionString = sprintf('sin(%s * n)', numToStr(omega0));
        case 3
            r = input('Enter r (real number, e.g., 0.9): ');
            validateNumeric(r);
            expressionString = sprintf('(%s).^abs(n)', numToStr(r));
        case 4
            k = input('Enter k (integer shift, e.g., 3): ');
            validateNumeric(k);
            k = round(k);
            expressionString = sprintf('d(n - (%s))', numToStr(k));
        case 5
            omega0 = input('Enter omega0 in radians (e.g., 0.1*pi): ');
            N = input('Enter N (positive integer length, e.g., 50): ');
            validateNumeric(omega0);
            validateNumeric(N);
            N = max(0, round(N));
            expressionString = sprintf('(cos(%s*n)).*(u(n) - u(n-(%s)))', numToStr(omega0), numToStr(N));
        case 6
            C = input('Enter constant C (e.g., 1): ');
            validateNumeric(C);
            expressionString = sprintf('%s + 0*n', numToStr(C));
        case 7
            A = input('Enter amplitude A (e.g., 2): ');
            r = input('Enter r (e.g., 0.8): ');
            validateNumeric(A);
            validateNumeric(r);
            expressionString = sprintf('(%s).*(%s).^n .* u(n)', numToStr(A), numToStr(r));
        case 8
            fprintf('Enter x[n] in terms of n. Helpers: u(n), d(n).\n');
            fprintf('Examples: (-0.5).^n.*u(n), sin(0.2*pi*n), (0.9).^abs(n), d(n-3)\n');
            expressionString = input('x[n] = ', 's');
            if isempty(strtrim(expressionString))
                error('No expression entered.');
            end
        otherwise
            error('Unexpected menu choice.');
    end
end

function validateNumeric(x)
    if isempty(x) || ~isscalar(x) || ~isfinite(x)
        error('Invalid numeric input.');
    end
end

function s = numToStr(x)
    % Robust numeric to string for constructing expressions
    if isinteger(x)
        s = num2str(double(x));
    else
        s = num2str(x, '%.15g');
    end
end
