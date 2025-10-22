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

%% Prompt for the expression
fprintf('Discrete-Time Signal Classifier (Energy / Power / Neither)\n');
fprintf('Enter x[n] in terms of n. Helpers: u(n), d(n).\n');
fprintf('Examples: (-0.5).^n.*u(n), sin(0.2*pi*n), (0.9).^abs(n), d(n-3)\n');
expressionString = input('x[n] = ', 's');

if isempty(strtrim(expressionString))
    error('No expression entered.');
end

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
