classdef SinogramViewer < handle
    % SinogramViewer - Specialized viewer for sinogram data
    %
    % This viewer provides tools for sinogram visualization, center of
    % rotation finding, and sinogram-specific analysis.
    %
    % Example:
    %   viewer = mufo.vis.SinogramViewer('Parent', panel);
    %   viewer.setData(sinogramData);
    %   center = viewer.findCenter();
    
    properties (Access = private)
        parent              % Parent container
        panel               % Main panel
        ax                  % Main axes
        imageHandle         % Image handle
        data                % Sinogram data
        currentSlice        % Current slice for 3D sinograms
        
        % UI Components
        controls            % Control panel
        sliceSlider         % Slice selector
        centerEdit          % Center of rotation edit
        angleInfo           % Angle information display
        
        % Analysis tools
        centerLine          % Center of rotation indicator
        sinoProfile         % Sinogram profile plot
        correlationPlot     % Correlation plot for center finding
        
        % Parameters
        angles              % Projection angles
        pixelSize           % Detector pixel size
        centerOfRotation    % Current center estimate
    end
    
    properties (Constant)
        DEFAULT_COLORMAP = 'gray';
    end
    
    events
        CenterChanged
        SliceChanged
    end
    
    methods
        function obj = SinogramViewer(varargin)
            % Constructor
            %
            % Parameters:
            %   'Parent' - Parent container
            %   'Position' - Position vector [x y w h]
            
            p = inputParser;
            addParameter(p, 'Parent', [], @ishandle);
            addParameter(p, 'Position', [0 0 1 1], @isnumeric);
            parse(p, varargin{:});
            
            obj.parent = p.Results.Parent;
            if isempty(obj.parent)
                obj.parent = figure('Name', 'Sinogram Viewer');
            end
            
            % Create UI
            obj.createUI(p.Results.Position);
            
            % Initialize
            obj.currentSlice = 1;
        end
        
        function setData(obj, data, angles)
            % Set sinogram data
            %
            % Inputs:
            %   data - Sinogram data (2D or 3D)
            %   angles - Projection angles in radians (optional)
            
            if ~isnumeric(data)
                error('SinogramViewer:InvalidData', ...
                    'Data must be numeric array');
            end
            
            obj.data = data;
            
            % Set angles
            if nargin < 3
                % Assume uniform distribution over 360 degrees
                nAngles = size(data, 1);
                obj.angles = linspace(0, 2*pi, nAngles+1);
                obj.angles(end) = [];
            else
                obj.angles = angles;
            end
            
            % Update UI
            obj.updateSliceControls();
            obj.updateDisplay();
            
            % Auto-detect center
            obj.findCenter('auto');
        end
        
        function center = findCenter(obj, method)
            % Find center of rotation
            %
            % Input:
            %   method - 'auto', 'correlation', 'entropy', 'manual'
            % Output:
            %   center - Center of rotation position
            
            if nargin < 2
                method = 'correlation';
            end
            
            sino = obj.getCurrentSinogram();
            
            switch lower(method)
                case 'auto'
                    % Try multiple methods and choose best
                    center1 = obj.findCenterCorrelation(sino);
                    center2 = obj.findCenterEntropy(sino);
                    
                    % Test both centers
                    score1 = obj.evaluateCenter(sino, center1);
                    score2 = obj.evaluateCenter(sino, center2);
                    
                    if score1 < score2
                        center = center1;
                    else
                        center = center2;
                    end
                    
                case 'correlation'
                    center = obj.findCenterCorrelation(sino);
                    
                case 'entropy'
                    center = obj.findCenterEntropy(sino);
                    
                case 'manual'
                    center = obj.manualCenterSelection();
                    
                otherwise
                    error('SinogramViewer:InvalidMethod', ...
                        'Unknown center finding method: %s', method);
            end
            
            obj.centerOfRotation = center;
            obj.updateCenterDisplay();
            
            notify(obj, 'CenterChanged');
        end
        
        function testReconstructionAtCenter(obj, center)
            % Show test reconstruction at given center
            %
            % Input:
            %   center - Center position to test
            
            if nargin < 2
                center = obj.centerOfRotation;
            end
            
            % Get current sinogram
            sino = obj.getCurrentSinogram();
            
            % Quick FBP reconstruction
            recon = obj.quickFBP(sino, center);
            
            % Show in new figure
            figure('Name', sprintf('Test Reconstruction - Center: %.2f', center));
            imagesc(recon);
            colormap gray;
            axis image;
            colorbar;
            title(sprintf('Center: %.2f', center));
        end
        
        function compareMultipleCenters(obj, centers)
            % Compare reconstructions at multiple centers
            %
            % Input:
            %   centers - Array of center positions to test
            
            if nargin < 2
                % Default: test around current center
                if isempty(obj.centerOfRotation)
                    obj.findCenter();
                end
                
                centers = obj.centerOfRotation + [-2, -1, 0, 1, 2];
            end
            
            sino = obj.getCurrentSinogram();
            
            % Create comparison figure
            figure('Name', 'Center of Rotation Comparison');
            n = length(centers);
            cols = ceil(sqrt(n));
            rows = ceil(n / cols);
            
            for i = 1:n
                subplot(rows, cols, i);
                
                % Quick reconstruction
                recon = obj.quickFBP(sino, centers(i));
                
                imagesc(recon);
                colormap gray;
                axis image off;
                title(sprintf('Center: %.1f', centers(i)));
            end
        end
        
        function analyzeSinogram(obj)
            % Comprehensive sinogram analysis
            
            sino = obj.getCurrentSinogram();
            
            % Create analysis figure
            fig = figure('Name', 'Sinogram Analysis', ...
                'Position', [100, 100, 1200, 800]);
            
            % Original sinogram
            subplot(2,3,1);
            imagesc(sino);
            colormap gray;
            xlabel('Detector Position');
            ylabel('Projection Angle');
            title('Original Sinogram');
            colorbar;
            
            % Angle-averaged profile
            subplot(2,3,2);
            avgProfile = mean(sino, 1);
            plot(avgProfile, 'b-', 'LineWidth', 2);
            xlabel('Detector Position');
            ylabel('Average Intensity');
            title('Detector Profile (Averaged)');
            grid on;
            
            % Add center line
            hold on;
            if ~isempty(obj.centerOfRotation)
                xline(obj.centerOfRotation, 'r--', 'LineWidth', 2);
                legend('Profile', 'Center', 'Location', 'best');
            end
            
            % Sinogram derivative
            subplot(2,3,3);
            sinoGrad = imgradient(sino);
            imagesc(sinoGrad);
            colormap hot;
            xlabel('Detector Position');
            ylabel('Projection Angle');
            title('Sinogram Gradient');
            colorbar;
            
            % FFT analysis
            subplot(2,3,4);
            fftSino = abs(fftshift(fft2(sino)));
            imagesc(log10(fftSino + 1));
            colormap jet;
            xlabel('Frequency X');
            ylabel('Frequency Y');
            title('2D FFT (log scale)');
            colorbar;
            
            % Projection consistency
            subplot(2,3,5);
            consistency = obj.checkProjectionConsistency(sino);
            plot(consistency, 'g-', 'LineWidth', 2);
            xlabel('Projection Number');
            ylabel('Consistency Score');
            title('Projection Consistency');
            grid on;
            
            % Statistics
            subplot(2,3,6);
            axis off;
            stats = obj.calculateSinogramStats(sino);
            
            statText = sprintf(['Sinogram Statistics:\n\n', ...
                'Size: %d x %d\n', ...
                'Min: %.3f\n', ...
                'Max: %.3f\n', ...
                'Mean: %.3f\n', ...
                'Std: %.3f\n', ...
                'Dynamic Range: %.1f dB\n', ...
                'Center: %.2f'], ...
                size(sino, 1), size(sino, 2), ...
                stats.min, stats.max, stats.mean, stats.std, ...
                stats.dynamicRange, obj.centerOfRotation);
            
            text(0.1, 0.9, statText, 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontName', 'FixedWidth');
        end
        
        function show(obj)
            % Show the viewer
            obj.panel.Visible = 'on';
        end
        
        function hide(obj)
            % Hide the viewer
            obj.panel.Visible = 'off';
        end
    end
    
    methods (Access = private)
        function createUI(obj, position)
            % Create UI components
            
            % Main panel
            obj.panel = uipanel(...
                'Parent', obj.parent, ...
                'Units', 'normalized', ...
                'Position', position, ...
                'BorderType', 'none');
            
            % Sinogram display (left)
            obj.ax = axes(...
                'Parent', obj.panel, ...
                'Units', 'normalized', ...
                'Position', [0.02, 0.1, 0.55, 0.85], ...
                'Box', 'on');
            
            % Control panel (right top)
            obj.controls = uipanel(...
                'Parent', obj.panel, ...
                'Title', 'Sinogram Controls', ...
                'Units', 'normalized', ...
                'Position', [0.58, 0.5, 0.4, 0.45]);
            
            % Analysis panel (right bottom)
            analysisPanel = uipanel(...
                'Parent', obj.panel, ...
                'Title', 'Center Finding', ...
                'Units', 'normalized', ...
                'Position', [0.58, 0.05, 0.4, 0.44]);
            
            % Create controls
            obj.createSinogramControls();
            obj.createCenterControls(analysisPanel);
            
            % Info display
            obj.createInfoDisplay();
        end
        
        function createSinogramControls(obj)
            % Create sinogram-specific controls
            
            yPos = 0.85;
            
            % Slice control (for 3D)
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Slice:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.2, 0.08], ...
                'HorizontalAlignment', 'left');
            
            obj.sliceSlider = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.25, yPos, 0.5, 0.08], ...
                'Min', 1, 'Max', 1, 'Value', 1, ...
                'Callback', @(src,~) obj.onSliceChange(src));
            
            % Slice number display
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', '1/1', ...
                'Units', 'normalized', ...
                'Position', [0.77, yPos, 0.2, 0.08], ...
                'Tag', 'SliceText');
            
            % Display options
            yPos = yPos - 0.15;
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'checkbox', ...
                'String', 'Show Center Line', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.45, 0.08], ...
                'Value', 1, ...
                'Callback', @(src,~) obj.toggleCenterLine(src));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'checkbox', ...
                'String', 'Log Scale', ...
                'Units', 'normalized', ...
                'Position', [0.5, yPos, 0.45, 0.08], ...
                'Callback', @(src,~) obj.toggleLogScale(src));
            
            % Analysis buttons
            yPos = yPos - 0.15;
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Analyze', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.42, 0.1], ...
                'Callback', @(~,~) obj.analyzeSinogram());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Compare Centers', ...
                'Units', 'normalized', ...
                'Position', [0.53, yPos, 0.42, 0.1], ...
                'Callback', @(~,~) obj.compareMultipleCenters());
            
            % View options
            yPos = yPos - 0.15;
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Show Profile', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.42, 0.1], ...
                'Callback', @(~,~) obj.showSinogramProfile());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Show Spectrum', ...
                'Units', 'normalized', ...
                'Position', [0.53, yPos, 0.42, 0.1], ...
                'Callback', @(~,~) obj.showFrequencySpectrum());
        end
        
        function createCenterControls(obj, parent)
            % Create center of rotation controls
            
            % Method selection
            uicontrol(...
                'Parent', parent, ...
                'Style', 'text', ...
                'String', 'Method:', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.85, 0.25, 0.08], ...
                'HorizontalAlignment', 'left');
            
            uicontrol(...
                'Parent', parent, ...
                'Style', 'popupmenu', ...
                'String', {'Auto', 'Correlation', 'Entropy', 'Manual'}, ...
                'Units', 'normalized', ...
                'Position', [0.3, 0.85, 0.65, 0.08], ...
                'Callback', @(src,~) obj.onMethodChange(src));
            
            % Current center
            uicontrol(...
                'Parent', parent, ...
                'Style', 'text', ...
                'String', 'Center:', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.7, 0.25, 0.08], ...
                'HorizontalAlignment', 'left');
            
            obj.centerEdit = uicontrol(...
                'Parent', parent, ...
                'Style', 'edit', ...
                'String', '', ...
                'Units', 'normalized', ...
                'Position', [0.3, 0.7, 0.3, 0.08], ...
                'Callback', @(src,~) obj.onCenterEdit(src));
            
            uicontrol(...
                'Parent', parent, ...
                'Style', 'pushbutton', ...
                'String', 'Test', ...
                'Units', 'normalized', ...
                'Position', [0.65, 0.7, 0.3, 0.08], ...
                'Callback', @(~,~) obj.testReconstructionAtCenter());
            
            % Fine adjustment
            uicontrol(...
                'Parent', parent, ...
                'Style', 'text', ...
                'String', 'Fine Adjust:', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.55, 0.9, 0.08], ...
                'HorizontalAlignment', 'left');
            
            % Adjustment buttons
            adjustments = {'-1', '-0.5', '-0.1', '+0.1', '+0.5', '+1'};
            btnWidth = 0.15;
            
            for i = 1:length(adjustments)
                xPos = 0.05 + (i-1) * btnWidth;
                uicontrol(...
                    'Parent', parent, ...
                    'Style', 'pushbutton', ...
                    'String', adjustments{i}, ...
                    'Units', 'normalized', ...
                    'Position', [xPos, 0.45, btnWidth-0.01, 0.08], ...
                    'Callback', @(~,~) obj.adjustCenter(str2double(adjustments{i})));
            end
            
            % Search range
            uicontrol(...
                'Parent', parent, ...
                'Style', 'text', ...
                'String', 'Search Range:', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.3, 0.9, 0.08], ...
                'HorizontalAlignment', 'left');
            
            uicontrol(...
                'Parent', parent, ...
                'Style', 'pushbutton', ...
                'String', 'Scan ±5 pixels', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.2, 0.42, 0.08], ...
                'Callback', @(~,~) obj.scanCenterRange(5));
            
            uicontrol(...
                'Parent', parent, ...
                'Style', 'pushbutton', ...
                'String', 'Scan ±10 pixels', ...
                'Units', 'normalized', ...
                'Position', [0.53, 0.2, 0.42, 0.08], ...
                'Callback', @(~,~) obj.scanCenterRange(10));
            
            % Optimization
            uicontrol(...
                'Parent', parent, ...
                'Style', 'pushbutton', ...
                'String', 'Optimize Center', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.05, 0.9, 0.1], ...
                'FontWeight', 'bold', ...
                'Callback', @(~,~) obj.optimizeCenter());
        end
        
        function createInfoDisplay(obj)
            % Create information display
            obj.angleInfo = uicontrol(...
                'Parent', obj.panel, ...
                'Style', 'text', ...
                'String', 'Angles: 0° - 360° | Projections: 0', ...
                'Units', 'normalized', ...
                'Position', [0.02, 0.02, 0.55, 0.06], ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor', [0.95, 0.95, 0.95], ...
                'FontName', 'FixedWidth');
        end
        
        function sino = getCurrentSinogram(obj)
            % Get current sinogram slice
            if ndims(obj.data) == 3
                sino = squeeze(obj.data(:, :, obj.currentSlice));
            else
                sino = obj.data;
            end
        end
        
        function updateDisplay(obj)
            % Update sinogram display
            if isempty(obj.data)
                return;
            end
            
            sino = obj.getCurrentSinogram();
            
            % Display sinogram
            if isempty(obj.imageHandle) || ~isvalid(obj.imageHandle)
                obj.imageHandle = imagesc(obj.ax, sino);
                colormap(obj.ax, obj.DEFAULT_COLORMAP);
                xlabel(obj.ax, 'Detector Position');
                ylabel(obj.ax, 'Projection Number');
                
                % Add center line
                hold(obj.ax, 'on');
                obj.centerLine = xline(obj.ax, 0, 'r--', 'LineWidth', 2);
                obj.centerLine.Visible = 'off';
                hold(obj.ax, 'off');
            else
                obj.imageHandle.CData = sino;
            end
            
            % Update title
            if ndims(obj.data) == 3
                title(obj.ax, sprintf('Sinogram - Slice %d/%d', ...
                    obj.currentSlice, size(obj.data, 3)));
            else
                title(obj.ax, 'Sinogram');
            end
            
            % Update info
            obj.updateInfo();
        end
        
        function updateInfo(obj)
            % Update information display
            if isempty(obj.data)
                return;
            end
            
            [nAngles, nDetectors] = size(obj.getCurrentSinogram());
            
            if ~isempty(obj.angles)
                angleRange = sprintf('%.1f° - %.1f°', ...
                    rad2deg(min(obj.angles)), rad2deg(max(obj.angles)));
            else
                angleRange = '0° - 360°';
            end
            
            infoStr = sprintf('Angles: %s | Projections: %d | Detectors: %d', ...
                angleRange, nAngles, nDetectors);
            
            if ~isempty(obj.centerOfRotation)
                infoStr = sprintf('%s | Center: %.2f', ...
                    infoStr, obj.centerOfRotation);
            end
            
            obj.angleInfo.String = infoStr;
        end
        
        function updateSliceControls(obj)
            % Update slice controls for 3D data
            if ndims(obj.data) == 3
                nSlices = size(obj.data, 3);
                
                obj.sliceSlider.Max = nSlices;
                obj.sliceSlider.Value = obj.currentSlice;
                obj.sliceSlider.Enable = 'on';
                
                if nSlices > 1
                    obj.sliceSlider.SliderStep = [1/(nSlices-1), 10/(nSlices-1)];
                end
                
                % Update text
                sliceText = findobj(obj.controls, 'Tag', 'SliceText');
                if ~isempty(sliceText)
                    sliceText.String = sprintf('%d/%d', obj.currentSlice, nSlices);
                end
            else
                obj.sliceSlider.Enable = 'off';
            end
        end
        
        function updateCenterDisplay(obj)
            % Update center of rotation display
            if ~isempty(obj.centerOfRotation)
                obj.centerEdit.String = sprintf('%.2f', obj.centerOfRotation);
                
                % Update center line
                if ~isempty(obj.centerLine) && isvalid(obj.centerLine)
                    obj.centerLine.Value = obj.centerOfRotation;
                    if obj.centerLine.Visible == "on"
                        obj.centerLine.Visible = 'on';
                    end
                end
            end
            
            obj.updateInfo();
        end
        
        function center = findCenterCorrelation(obj, sino)
            % Find center using correlation method
            %
            % Based on correlation between projections at 0° and 180°
            
            [nAngles, nDetectors] = size(sino);
            
            % Find opposite projections
            if nAngles >= 180
                % Find projections closest to 0° and 180°
                idx1 = 1;
                idx2 = round(nAngles / 2);
            else
                warning('SinogramViewer:InsufficientAngles', ...
                    'Not enough angles for correlation method');
                center = nDetectors / 2;
                return;
            end
            
            % Get projections
            proj1 = sino(idx1, :);
            proj2 = sino(idx2, :);
            
            % Flip second projection
            proj2_flipped = fliplr(proj2);
            
            % Cross-correlation
            correlation = xcorr(proj1, proj2_flipped);
            [~, maxIdx] = max(correlation);
            
            % Calculate center
            shift = maxIdx - length(proj1);
            center = nDetectors/2 + shift/2;
            
            % Refine with sub-pixel accuracy
            if maxIdx > 1 && maxIdx < length(correlation)
                % Parabolic fit around maximum
                y1 = correlation(maxIdx-1);
                y2 = correlation(maxIdx);
                y3 = correlation(maxIdx+1);
                
                delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
                center = center + delta/2;
            end
        end
        
        function center = findCenterEntropy(obj, sino)
            % Find center using entropy minimization
            
            [~, nDetectors] = size(sino);
            
            % Search range
            searchRange = round(nDetectors * 0.1); % 10% of detector width
            centerGuess = nDetectors / 2;
            
            centers = centerGuess + (-searchRange:0.5:searchRange);
            entropies = zeros(size(centers));
            
            % Test each center
            for i = 1:length(centers)
                % Quick reconstruction
                recon = obj.quickFBP(sino, centers(i));
                
                % Calculate image entropy
                recon_norm = mat2gray(recon);
                entropies(i) = entropy(recon_norm);
            end
            
            % Find minimum entropy
            [~, minIdx] = min(entropies);
            center = centers(minIdx);
            
            % Refine with finer search
            if minIdx > 1 && minIdx < length(centers)
                % Fine search around minimum
                fineRange = centers(minIdx+1) - centers(minIdx);
                fineCenters = center + (-fineRange:0.1:fineRange);
                fineEntropies = zeros(size(fineCenters));
                
                for i = 1:length(fineCenters)
                    recon = obj.quickFBP(sino, fineCenters(i));
                    recon_norm = mat2gray(recon);
                    fineEntropies(i) = entropy(recon_norm);
                end
                
                [~, minIdx] = min(fineEntropies);
                center = fineCenters(minIdx);
            end
        end
        
        function score = evaluateCenter(obj, sino, center)
            % Evaluate quality of center position
            
            % Quick reconstruction
            recon = obj.quickFBP(sino, center);
            
            % Multiple metrics
            metrics = struct();
            
            % 1. Image entropy (lower is better)
            recon_norm = mat2gray(recon);
            metrics.entropy = entropy(recon_norm);
            
            % 2. Edge strength (higher is better)
            [Gmag, ~] = imgradient(recon);
            metrics.edgeStrength = mean(Gmag(:));
            
            % 3. Symmetry (lower is better)
            recon_left = recon(:, 1:floor(end/2));
            recon_right = fliplr(recon(:, ceil(end/2)+1:end));
            minSize = min(size(recon_left, 2), size(recon_right, 2));
            metrics.asymmetry = mean(abs(recon_left(:, 1:minSize) - ...
                                        recon_right(:, 1:minSize)), 'all');
            
            % Combined score (lower is better)
            score = metrics.entropy + metrics.asymmetry / metrics.edgeStrength;
        end
        
        function recon = quickFBP(obj, sino, center)
            % Quick filtered backprojection for center testing
            
            [nAngles, nDetectors] = size(sino);
            
            % Create filter (Ram-Lak)
            filter = obj.createRamLakFilter(nDetectors);
            
            % Apply filter to each projection
            sinoFiltered = real(ifft(fft(sino, [], 2) .* ...
                                    repmat(filter, nAngles, 1), [], 2));
            
            % Backprojection
            outputSize = nDetectors;
            recon = zeros(outputSize, outputSize);
            
            [X, Y] = meshgrid(1:outputSize, 1:outputSize);
            xCenter = outputSize/2;
            yCenter = outputSize/2;
            
            for i = 1:nAngles
                angle = obj.angles(i);
                
                % Projection coordinates
                t = (X - xCenter) * cos(angle) + (Y - yCenter) * sin(angle);
                tCoord = t + center;
                
                % Interpolate
                proj = sinoFiltered(i, :);
                projInterp = interp1(1:nDetectors, proj, tCoord, 'linear', 0);
                
                % Accumulate
                recon = recon + projInterp;
            end
            
            recon = recon * pi / nAngles;
        end
        
        function filter = createRamLakFilter(obj, n)
            % Create Ram-Lak filter for FBP
            
            % Frequency axis
            freq = [0:n/2-1, -n/2:-1] / n;
            
            % Ram-Lak filter
            filter = abs(freq);
            
            % Apply Hamming window to reduce ringing
            hamming = 0.54 + 0.46 * cos(2*pi*freq);
            filter = filter .* hamming;
            
            filter = fftshift(filter);
        end
        
        function center = manualCenterSelection(obj)
            % Manual center selection with mouse
            
            title(obj.ax, 'Click to select center of rotation');
            
            % Get point
            [x, ~] = ginput(1);
            center = x;
            
            title(obj.ax, '');
        end
        
        function consistency = checkProjectionConsistency(obj, sino)
            % Check consistency between adjacent projections
            
            nAngles = size(sino, 1);
            consistency = zeros(nAngles-1, 1);
            
            for i = 1:nAngles-1
                % Compare adjacent projections
                diff = abs(sino(i, :) - sino(i+1, :));
                consistency(i) = mean(diff);
            end
            
            % Normalize
            consistency = consistency / mean(consistency);
        end
        
        function stats = calculateSinogramStats(obj, sino)
            % Calculate sinogram statistics
            
            stats = struct();
            stats.min = min(sino(:));
            stats.max = max(sino(:));
            stats.mean = mean(sino(:));
            stats.std = std(sino(:));
            
            % Dynamic range in dB
            if stats.max > 0 && stats.min > 0
                stats.dynamicRange = 20 * log10(stats.max / stats.min);
            else
                stats.dynamicRange = inf;
            end
        end
        
        function showSinogramProfile(obj)
            % Show sinogram profiles
            
            sino = obj.getCurrentSinogram();
            
            figure('Name', 'Sinogram Profiles');
            
            % Horizontal profile (through angles)
            subplot(2,1,1);
            midDetector = round(size(sino, 2) / 2);
            profile = sino(:, midDetector);
            plot(rad2deg(obj.angles), profile, 'b-', 'LineWidth', 2);
            xlabel('Angle (degrees)');
            ylabel('Intensity');
            title(sprintf('Angular Profile at Detector %d', midDetector));
            grid on;
            
            % Vertical profile (through detectors)
            subplot(2,1,2);
            midAngle = round(size(sino, 1) / 2);
            profile = sino(midAngle, :);
            plot(profile, 'r-', 'LineWidth', 2);
            xlabel('Detector Position');
            ylabel('Intensity');
            title(sprintf('Detector Profile at Angle %.1f°', ...
                rad2deg(obj.angles(midAngle))));
            grid on;
            
            % Add center line
            if ~isempty(obj.centerOfRotation)
                hold on;
                xline(obj.centerOfRotation, 'k--', 'LineWidth', 2);
                legend('Profile', 'Center', 'Location', 'best');
            end
        end
        
        function showFrequencySpectrum(obj)
            % Show frequency spectrum analysis
            
            sino = obj.getCurrentSinogram();
            
            figure('Name', 'Sinogram Frequency Analysis');
            
            % 2D FFT
            subplot(2,2,1);
            fftSino = fftshift(fft2(sino));
            imagesc(log10(abs(fftSino) + 1));
            colormap jet;
            colorbar;
            xlabel('Frequency X');
            ylabel('Frequency Y');
            title('2D FFT Magnitude (log)');
            
            % Radial average
            subplot(2,2,2);
            [radialAvg, radialFreq] = obj.radialAverage(abs(fftSino));
            loglog(radialFreq, radialAvg, 'b-', 'LineWidth', 2);
            xlabel('Radial Frequency');
            ylabel('Magnitude');
            title('Radial Power Spectrum');
            grid on;
            
            % Horizontal frequency profile
            subplot(2,2,3);
            centerY = round(size(fftSino, 1) / 2);
            hProfile = abs(fftSino(centerY, :));
            plot(hProfile, 'g-', 'LineWidth', 2);
            xlabel('Horizontal Frequency');
            ylabel('Magnitude');
            title('Horizontal Frequency Profile');
            grid on;
            
            % Vertical frequency profile
            subplot(2,2,4);
            centerX = round(size(fftSino, 2) / 2);
            vProfile = abs(fftSino(:, centerX));
            plot(vProfile, 'r-', 'LineWidth', 2);
            xlabel('Vertical Frequency');
            ylabel('Magnitude');
            title('Vertical Frequency Profile');
            grid on;
        end
        
        function [avg, freq] = radialAverage(obj, img)
            % Compute radial average of 2D image
            
            [ny, nx] = size(img);
            centerX = (nx + 1) / 2;
            centerY = (ny + 1) / 2;
            
            [X, Y] = meshgrid(1:nx, 1:ny);
            R = sqrt((X - centerX).^2 + (Y - centerY).^2);
            
            maxR = min(centerX, centerY);
            nBins = round(maxR);
            
            avg = zeros(nBins, 1);
            freq = (0:nBins-1)';
            
            for i = 1:nBins
                mask = (R >= i-1) & (R < i);
                avg(i) = mean(img(mask));
            end
        end
        
        function scanCenterRange(obj, range)
            % Scan a range of center positions
            
            sino = obj.getCurrentSinogram();
            
            if isempty(obj.centerOfRotation)
                obj.findCenter();
            end
            
            centers = obj.centerOfRotation + (-range:0.5:range);
            scores = zeros(size(centers));
            
            % Progress bar
            h = waitbar(0, 'Scanning center positions...');
            
            for i = 1:length(centers)
                scores(i) = obj.evaluateCenter(sino, centers(i));
                waitbar(i/length(centers), h);
            end
            
            close(h);
            
            % Find best center
            [~, minIdx] = min(scores);
            bestCenter = centers(minIdx);
            
            % Show results
            figure('Name', 'Center Scan Results');
            plot(centers, scores, 'b-', 'LineWidth', 2);
            hold on;
            plot(bestCenter, scores(minIdx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            xlabel('Center Position');
            ylabel('Quality Score (lower is better)');
            title('Center of Rotation Scan');
            grid on;
            legend('Score', 'Best Center', 'Location', 'best');
            
            % Update center
            obj.centerOfRotation = bestCenter;
            obj.updateCenterDisplay();
        end
        
        function optimizeCenter(obj)
            % Optimize center using multiple methods
            
            sino = obj.getCurrentSinogram();
            
            % Try different methods
            methods = {'correlation', 'entropy'};
            centers = zeros(length(methods), 1);
            
            h = waitbar(0, 'Optimizing center...');
            
            for i = 1:length(methods)
                centers(i) = obj.findCenter(methods{i});
                waitbar(i/length(methods), h);
            end
            
            close(h);
            
            % Evaluate all centers
            scores = zeros(size(centers));
            for i = 1:length(centers)
                scores(i) = obj.evaluateCenter(sino, centers(i));
            end
            
            % Choose best
            [~, bestIdx] = min(scores);
            obj.centerOfRotation = centers(bestIdx);
            obj.updateCenterDisplay();
            
            % Show comparison
            obj.compareMultipleCenters(centers);
        end
        
        % Callback functions
        function onSliceChange(obj, src)
            % Handle slice change
            obj.currentSlice = round(src.Value);
            obj.updateDisplay();
            obj.updateSliceControls();
            notify(obj, 'SliceChanged');
        end
        
        function onMethodChange(obj, src)
            % Handle center finding method change
            methods = {'auto', 'correlation', 'entropy', 'manual'};
            method = methods{src.Value};
            obj.findCenter(method);
        end
        
        function onCenterEdit(obj, src)
            % Handle manual center input
            val = str2double(src.String);
            if ~isnan(val) && val > 0 && val < size(obj.data, 2)
                obj.centerOfRotation = val;
                obj.updateCenterDisplay();
                notify(obj, 'CenterChanged');
            else
                % Reset to current value
                src.String = sprintf('%.2f', obj.centerOfRotation);
            end
        end
        
        function adjustCenter(obj, delta)
            % Adjust center by delta
            if ~isempty(obj.centerOfRotation)
                obj.centerOfRotation = obj.centerOfRotation + delta;
                obj.updateCenterDisplay();
                notify(obj, 'CenterChanged');
            end
        end
        
        function toggleCenterLine(obj, src)
            % Toggle center line visibility
            if ~isempty(obj.centerLine) && isvalid(obj.centerLine)
                if src.Value
                    obj.centerLine.Visible = 'on';
                else
                    obj.centerLine.Visible = 'off';
                end
            end
        end
        
        function toggleLogScale(obj, src)
            % Toggle logarithmic scale
            if src.Value
                sino = obj.getCurrentSinogram();
                sino_log = log10(sino + 1);
                obj.imageHandle.CData = sino_log;
            else
                obj.updateDisplay();
            end
        end
    end
end