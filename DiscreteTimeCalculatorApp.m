classdef DiscreteTimeCalculatorApp < matlab.apps.AppBase
    % Discrete-time LCCDE Calculator (uifigure-based)
    % - Direct method: zero-input (homogeneous), zero-state (particular), total
    % - Indirect method: unilateral z-transform via impulse response + convolution

    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        RootGrid                   matlab.ui.container.GridLayout

        % Left controls
        ControlGrid                matlab.ui.container.GridLayout
        ALabel                     matlab.ui.control.Label
        AEditField                 matlab.ui.control.EditField
        BLabel                     matlab.ui.control.Label
        BEditField                 matlab.ui.control.EditField
        NLabel                     matlab.ui.control.Label
        NEditField                 matlab.ui.control.NumericEditField
        MethodLabel                matlab.ui.control.Label
        MethodDropDown             matlab.ui.control.DropDown

        InputTypeLabel             matlab.ui.control.Label
        InputTypeDropDown          matlab.ui.control.DropDown
        XVectorLabel               matlab.ui.control.Label
        XVectorEditField           matlab.ui.control.EditField
        AmpLabel                   matlab.ui.control.Label
        AmpEditField               matlab.ui.control.NumericEditField
        RLabel                     matlab.ui.control.Label
        REditField                 matlab.ui.control.NumericEditField
        OmegaLabel                 matlab.ui.control.Label
        OmegaEditField             matlab.ui.control.NumericEditField
        PhiLabel                   matlab.ui.control.Label
        PhiEditField               matlab.ui.control.NumericEditField

        ICTypeLabel                matlab.ui.control.Label
        ICTypeDropDown             matlab.ui.control.DropDown
        ICLabel                    matlab.ui.control.Label
        ICEditField                matlab.ui.control.EditField

        ShowComponentsCheckBox     matlab.ui.control.CheckBox
        ComputeButton              matlab.ui.control.Button
        ExportButton               matlab.ui.control.Button
        StatusLabel                matlab.ui.control.Label

        % Right plots
        TabGroup                   matlab.ui.container.TabGroup
        XTab                       matlab.ui.container.Tab
        XAxes                      matlab.ui.control.UIAxes
        YTab                       matlab.ui.container.Tab
        YAxes                      matlab.ui.control.UIAxes
    end

    properties (Access = private)
        lastX double = []
        lastYZeroInput double = []
        lastYZeroState double = []
        lastYTotal double = []
        lastImpulseResponse double = []
    end

    methods (Access = private)
        function onInputTypeChanged(app, ~)
            t = app.InputTypeDropDown.Value;
            % Hide all optional fields by default
            app.XVectorLabel.Visible = false; app.XVectorEditField.Visible = false;
            app.AmpLabel.Visible = false; app.AmpEditField.Visible = false;
            app.RLabel.Visible = false; app.REditField.Visible = false;
            app.OmegaLabel.Visible = false; app.OmegaEditField.Visible = false;
            app.PhiLabel.Visible = false; app.PhiEditField.Visible = false;

            switch t
                case 'Custom vector'
                    app.XVectorLabel.Visible = true; app.XVectorEditField.Visible = true;
                case 'Impulse'
                    app.AmpLabel.Visible = true; app.AmpEditField.Visible = true;
                case 'Step'
                    app.AmpLabel.Visible = true; app.AmpEditField.Visible = true;
                case 'Exponential A*r^n'
                    app.AmpLabel.Visible = true; app.AmpEditField.Visible = true;
                    app.RLabel.Visible = true; app.REditField.Visible = true;
                case 'Sine A*sin(ωn+φ)'
                    app.AmpLabel.Visible = true; app.AmpEditField.Visible = true;
                    app.OmegaLabel.Visible = true; app.OmegaEditField.Visible = true;
                    app.PhiLabel.Visible = true; app.PhiEditField.Visible = true;
            end
        end

        function onCompute(app, ~)
            try
                app.showStatus('Computing...', false);
                drawnow;

                [a, b, N] = app.parseCoreParams();

                % Build input sequence
                inputType = app.InputTypeDropDown.Value;
                A = app.AmpEditField.Value;
                r = app.REditField.Value;
                omega = app.OmegaEditField.Value;
                phi = app.PhiEditField.Value;
                x = app.buildInputVector(inputType, N, A, r, omega, phi, app.XVectorEditField.Value);

                % Initial conditions
                icMode = app.ICTypeDropDown.Value; % 'All zero' | 'y[-1..-N]' | 'y[0..N-1]'
                orderA = max(length(a) - 1, 0);
                icVector = zeros(1, orderA);
                initialMode = 'left';
                if orderA > 0
                    switch icMode
                        case 'All zero'
                            icVector = zeros(1, orderA);
                            initialMode = 'left';
                        case 'y[-1..-N]'
                            icVector = app.parseVector(app.ICEditField.Value, 'Initial conditions y[-1..-N]');
                            icVector = app.ensureLength(icVector, orderA);
                            initialMode = 'left';
                        case 'y[0..N-1]'
                            icVector = app.parseVector(app.ICEditField.Value, 'Initial conditions y[0..N-1]');
                            icVector = app.ensureLength(icVector, orderA);
                            initialMode = 'right';
                    end
                end

                method = app.MethodDropDown.Value; % 'Direct' | 'Indirect (z-transform)'

                % Direct method zero-state
                yZS_direct = app.zeroStateResponse(a, b, x, N);
                % Zero-input from initial conditions
                yZI = app.zeroInputResponse(a, icVector, N, initialMode);

                if strcmp(method, 'Direct')
                    yZS = yZS_direct;
                else
                    % Indirect: compute impulse response and convolve
                    h = app.impulseResponse(a, b, N);
                    yZS = conv(x, h);
                    yZS = yZS(1:N);
                    app.lastImpulseResponse = h;
                end

                yTOTAL = yZS + yZI;

                % Store last results
                app.lastX = x; app.lastYZeroState = yZS; app.lastYZeroInput = yZI; app.lastYTotal = yTOTAL;

                % Plots
                n = 0:(N-1);
                % x[n]
                cla(app.XAxes);
                stem(app.XAxes, n, x);
                grid(app.XAxes, 'on'); xlabel(app.XAxes, 'n'); ylabel(app.XAxes, 'x[n]');
                title(app.XAxes, 'Input sequence x[n]');

                % y[n]
                cla(app.YAxes);
                hold(app.YAxes, 'on');
                if app.ShowComponentsCheckBox.Value
                    stem(app.YAxes, n, yZI, 'DisplayName', 'Zero-input (homogeneous)');
                    stem(app.YAxes, n, yZS, 'DisplayName', 'Zero-state (particular)');
                end
                stem(app.YAxes, n, yTOTAL, 'DisplayName', 'Total y[n]');
                hold(app.YAxes, 'off');
                grid(app.YAxes, 'on'); xlabel(app.YAxes, 'n'); ylabel(app.YAxes, 'y[n]');
                title(app.YAxes, sprintf('Output (%s)', method));
                legend(app.YAxes, 'Location', 'best');

                app.showStatus('Done.', false);
            catch ME
                app.showStatus(ME.message, true);
            end
        end

        function onExport(app, ~)
            try
                if isempty(app.lastX) || isempty(app.lastYTotal)
                    error('Nothing to export. Compute first.');
                end
                N = numel(app.lastX);
                n = (0:N-1).';
                x = app.lastX(:);
                yzi = app.lastYZeroInput(:);
                yzs = app.lastYZeroState(:);
                y = app.lastYTotal(:);
                h = app.lastImpulseResponse(:);
                if isempty(h)
                    h = nan(N,1);
                else
                    if numel(h) < N
                        h = [h(:); nan(N - numel(h), 1)];
                    elseif numel(h) > N
                        h = h(1:N);
                    end
                end

                data = [n, x, yzs, yzi, y, h];
                [file, path] = uiputfile({'*.csv', 'CSV file (*.csv)'}, 'Export sequences as CSV', 'dtime_calc_export.csv');
                if isequal(file, 0)
                    return;
                end
                outFile = fullfile(path, file);
                % Prefer writetable with headers for portability
                T = table(n, x, yzs, yzi, y, h, 'VariableNames', ...
                    {'n','x','y_zero_state','y_zero_input','y_total','h'});
                writetable(T, outFile);
                app.showStatus(sprintf('Exported: %s', outFile), false);
            catch ME
                app.showStatus(ME.message, true);
            end
        end

        function [a, b, N] = parseCoreParams(app)
            a = app.parseVector(app.AEditField.Value, 'Denominator a');
            b = app.parseVector(app.BEditField.Value, 'Numerator b');
            if isempty(a)
                error('Denominator a must not be empty.');
            end
            if a(1) == 0
                error('a(1) must be nonzero.');
            end
            if isempty(b)
                b = 0; % treat empty numerator as zero input feedthrough
            end
            N = max(1, floor(app.NEditField.Value));
            app.NEditField.Value = N;
        end

        function v = parseVector(~, s, name)
            if isstring(s); s = char(s); end
            if isempty(strtrim(s))
                v = [];
                return;
            end
            v = str2num(s); %#ok<ST2NM>
            if isempty(v)
                error('%s: invalid vector. Use MATLAB vector syntax, e.g., [1 -0.5 0.25].', name);
            end
            v = v(:).';
        end

        function v = ensureLength(~, v, L)
            if numel(v) < L
                v = [v(:).', zeros(1, L - numel(v))];
            elseif numel(v) > L
                v = v(1:L);
            end
        end

        function x = buildInputVector(app, type, N, A, r, omega, phi, xText)
            n = 0:(N-1);
            switch type
                case 'Custom vector'
                    x = app.parseVector(xText, 'x[n]');
                    if isempty(x)
                        error('x[n]: provide a non-empty vector for Custom vector input.');
                    end
                    if numel(x) < N
                        x = [x(:).', zeros(1, N - numel(x))];
                    elseif numel(x) > N
                        x = x(1:N);
                    else
                        x = x(:).';
                    end
                case 'Impulse'
                    x = zeros(1, N); x(1) = A;
                case 'Step'
                    x = A * ones(1, N);
                case 'Exponential A*r^n'
                    x = A * (r .^ n);
                case 'Sine A*sin(ωn+φ)'
                    x = A * sin(omega * n + phi);
                otherwise
                    error('Unknown input type.');
            end
        end

        function y = zeroStateResponse(~, a, b, x, N)
            a0 = a(1);
            orderA = max(length(a) - 1, 0);
            orderB = max(length(b) - 1, 0);
            y = zeros(1, N);
            for n = 0:(N-1)
                sumB = 0.0;
                for k = 0:orderB
                    if (n - k) >= 0
                        sumB = sumB + b(k+1) * x(n - k + 1);
                    end
                end
                sumA = 0.0;
                for k = 1:orderA
                    if (n - k) >= 0
                        sumA = sumA + a(k+1) * y(n - k + 1);
                    end
                end
                y(n+1) = (sumB - sumA) / a0;
            end
        end

        function y = zeroInputResponse(~, a, icVector, N, initialMode)
            a0 = a(1);
            orderA = max(length(a) - 1, 0);
            y = zeros(1, N);
            if orderA == 0
                return;
            end
            switch initialMode
                case 'left' % icVector corresponds to y(-1), y(-2), ..., y(-orderA)
                    for n = 0:(N-1)
                        sumA = 0.0;
                        for k = 1:orderA
                            if (n - k) >= 0
                                yprev = y(n - k + 1);
                            else
                                yprev = icVector(k); % y(-k)
                            end
                            sumA = sumA + a(k+1) * yprev;
                        end
                        y(n+1) = (-sumA) / a0;
                    end
                case 'right' % icVector corresponds to y(0), y(1), ..., y(orderA-1)
                    L = min(orderA, N);
                    if L > 0
                        y(1:L) = icVector(1:L);
                    end
                    for n = orderA:(N-1)
                        sumA = 0.0;
                        for k = 1:orderA
                            sumA = sumA + a(k+1) * y(n - k + 1);
                        end
                        y(n+1) = (-sumA) / a0;
                    end
                otherwise
                    error('Unknown initialMode: %s', initialMode);
            end
        end

        function h = impulseResponse(app, a, b, N)
            ximp = zeros(1, N); ximp(1) = 1;
            h = app.zeroStateResponse(a, b, ximp, N);
        end

        function showStatus(app, msg, isError)
            app.StatusLabel.Text = char(string(msg));
            if isError
                app.StatusLabel.FontColor = [0.8 0 0];
            else
                app.StatusLabel.FontColor = [0.1 0.5 0.1];
            end
        end
    end

    methods (Access = private)
        function createComponents(app)
            % UIFigure and layout
            app.UIFigure = uifigure('Name', 'Discrete Time Calculator');
            app.UIFigure.Position(3:4) = [1100, 650];
            app.UIFigure.CloseRequestFcn = @(~, ~) delete(app);

            app.RootGrid = uigridlayout(app.UIFigure, [1, 2]);
            app.RootGrid.ColumnWidth = {380, '1x'};

            % Left controls grid
            app.ControlGrid = uigridlayout(app.RootGrid, [17, 2]);
            app.ControlGrid.Layout.Row = 1; app.ControlGrid.Layout.Column = 1;
            app.ControlGrid.RowHeight = repmat({28}, 1, 17);
            app.ControlGrid.ColumnWidth = {'fit','1x'};

            % Coefficients
            app.ALabel = uilabel(app.ControlGrid); app.ALabel.Text = 'Denominator a (e.g., [1 -0.5 0.25])';
            app.ALabel.Layout.Row = 1; app.ALabel.Layout.Column = [1 2];
            app.AEditField = uieditfield(app.ControlGrid, 'text');
            app.AEditField.Layout.Row = 2; app.AEditField.Layout.Column = [1 2];
            app.AEditField.Value = '[1 -0.5 0.25]';

            app.BLabel = uilabel(app.ControlGrid); app.BLabel.Text = 'Numerator b (e.g., [0.2 0.1])';
            app.BLabel.Layout.Row = 3; app.BLabel.Layout.Column = [1 2];
            app.BEditField = uieditfield(app.ControlGrid, 'text');
            app.BEditField.Layout.Row = 4; app.BEditField.Layout.Column = [1 2];
            app.BEditField.Value = '[0.2 0.1]';

            app.NLabel = uilabel(app.ControlGrid); app.NLabel.Text = 'Simulation length N';
            app.NLabel.Layout.Row = 5; app.NLabel.Layout.Column = 1;
            app.NEditField = uieditfield(app.ControlGrid, 'numeric');
            app.NEditField.Layout.Row = 5; app.NEditField.Layout.Column = 2;
            app.NEditField.Limits = [1 Inf]; app.NEditField.RoundFractionalValues = 'on'; app.NEditField.Value = 64;

            app.MethodLabel = uilabel(app.ControlGrid); app.MethodLabel.Text = 'Method';
            app.MethodLabel.Layout.Row = 6; app.MethodLabel.Layout.Column = 1;
            app.MethodDropDown = uidropdown(app.ControlGrid, 'Items', {'Direct', 'Indirect (z-transform)'});
            app.MethodDropDown.Layout.Row = 6; app.MethodDropDown.Layout.Column = 2;

            % Input controls
            app.InputTypeLabel = uilabel(app.ControlGrid); app.InputTypeLabel.Text = 'Input x[n] type';
            app.InputTypeLabel.Layout.Row = 7; app.InputTypeLabel.Layout.Column = 1;
            app.InputTypeDropDown = uidropdown(app.ControlGrid, 'Items', { ...
                'Custom vector', 'Impulse', 'Step', 'Exponential A*r^n', 'Sine A*sin(ωn+φ)'});
            app.InputTypeDropDown.Layout.Row = 7; app.InputTypeDropDown.Layout.Column = 2;
            app.InputTypeDropDown.Value = 'Step';

            app.XVectorLabel = uilabel(app.ControlGrid); app.XVectorLabel.Text = 'x[n] (vector)';
            app.XVectorLabel.Layout.Row = 8; app.XVectorLabel.Layout.Column = 1;
            app.XVectorEditField = uieditfield(app.ControlGrid, 'text');
            app.XVectorEditField.Layout.Row = 8; app.XVectorEditField.Layout.Column = 2;
            app.XVectorEditField.Value = '[ ]';

            app.AmpLabel = uilabel(app.ControlGrid); app.AmpLabel.Text = 'Amplitude A';
            app.AmpLabel.Layout.Row = 9; app.AmpLabel.Layout.Column = 1;
            app.AmpEditField = uieditfield(app.ControlGrid, 'numeric');
            app.AmpEditField.Layout.Row = 9; app.AmpEditField.Layout.Column = 2; app.AmpEditField.Value = 1;

            app.RLabel = uilabel(app.ControlGrid); app.RLabel.Text = 'Base r (exp)';
            app.RLabel.Layout.Row = 10; app.RLabel.Layout.Column = 1;
            app.REditField = uieditfield(app.ControlGrid, 'numeric');
            app.REditField.Layout.Row = 10; app.REditField.Layout.Column = 2; app.REditField.Value = 0.9;

            app.OmegaLabel = uilabel(app.ControlGrid); app.OmegaLabel.Text = 'ω (rad/sample)';
            app.OmegaLabel.Layout.Row = 11; app.OmegaLabel.Layout.Column = 1;
            app.OmegaEditField = uieditfield(app.ControlGrid, 'numeric');
            app.OmegaEditField.Layout.Row = 11; app.OmegaEditField.Layout.Column = 2; app.OmegaEditField.Value = pi/6;

            app.PhiLabel = uilabel(app.ControlGrid); app.PhiLabel.Text = 'φ (rad)';
            app.PhiLabel.Layout.Row = 12; app.PhiLabel.Layout.Column = 1;
            app.PhiEditField = uieditfield(app.ControlGrid, 'numeric');
            app.PhiEditField.Layout.Row = 12; app.PhiEditField.Layout.Column = 2; app.PhiEditField.Value = 0;

            % Initial conditions
            app.ICTypeLabel = uilabel(app.ControlGrid); app.ICTypeLabel.Text = 'Initial conditions type';
            app.ICTypeLabel.Layout.Row = 13; app.ICTypeLabel.Layout.Column = 1;
            app.ICTypeDropDown = uidropdown(app.ControlGrid, 'Items', {'All zero', 'y[-1..-N]', 'y[0..N-1]'});
            app.ICTypeDropDown.Layout.Row = 13; app.ICTypeDropDown.Layout.Column = 2;

            app.ICLabel = uilabel(app.ControlGrid); app.ICLabel.Text = 'IC vector (match order)';
            app.ICLabel.Layout.Row = 14; app.ICLabel.Layout.Column = 1;
            app.ICEditField = uieditfield(app.ControlGrid, 'text');
            app.ICEditField.Layout.Row = 14; app.ICEditField.Layout.Column = 2;
            app.ICEditField.Value = '';

            app.ShowComponentsCheckBox = uicheckbox(app.ControlGrid, 'Text', 'Show homogeneous/particular');
            app.ShowComponentsCheckBox.Layout.Row = 15; app.ShowComponentsCheckBox.Layout.Column = [1 2];
            app.ShowComponentsCheckBox.Value = true;

            btnGrid = uigridlayout(app.ControlGrid, [1, 2]);
            btnGrid.Layout.Row = 16; btnGrid.Layout.Column = [1 2];
            btnGrid.ColumnWidth = {'1x','1x'};
            app.ComputeButton = uibutton(btnGrid, 'push', 'Text', 'Compute');
            app.ComputeButton.Layout.Row = 1; app.ComputeButton.Layout.Column = 1;
            app.ExportButton = uibutton(btnGrid, 'push', 'Text', 'Export CSV');
            app.ExportButton.Layout.Row = 1; app.ExportButton.Layout.Column = 2;

            app.StatusLabel = uilabel(app.ControlGrid);
            app.StatusLabel.Layout.Row = 17; app.StatusLabel.Layout.Column = [1 2];
            app.StatusLabel.HorizontalAlignment = 'center';
            app.StatusLabel.Text = '';

            % Right tabs and axes
            app.TabGroup = uitabgroup(app.RootGrid);
            app.TabGroup.Layout.Row = 1; app.TabGroup.Layout.Column = 2;

            app.XTab = uitab(app.TabGroup, 'Title', 'x[n]');
            app.XAxes = uiaxes(app.XTab);
            app.XAxes.Position = [30 30 650 520];

            app.YTab = uitab(app.TabGroup, 'Title', 'y[n]');
            app.YAxes = uiaxes(app.YTab);
            app.YAxes.Position = [30 30 650 520];

            % Wiring callbacks
            app.InputTypeDropDown.ValueChangedFcn = @(src, evt) app.onInputTypeChanged(src);
            app.ComputeButton.ButtonPushedFcn = @(src, evt) app.onCompute(src);
            app.ExportButton.ButtonPushedFcn = @(src, evt) app.onExport(src);

            % Initialize visibility
            app.onInputTypeChanged();

            % Initial placeholder plots
            cla(app.XAxes); cla(app.YAxes);
            title(app.XAxes, 'Input sequence x[n]'); xlabel(app.XAxes, 'n'); ylabel(app.XAxes, 'x[n]'); grid(app.XAxes, 'on');
            title(app.YAxes, 'Output y[n]'); xlabel(app.YAxes, 'n'); ylabel(app.YAxes, 'y[n]'); grid(app.YAxes, 'on');
        end
    end

    methods (Access = public)
        function app = DiscreteTimeCalculatorApp
            createComponents(app);
        end
    end
end
