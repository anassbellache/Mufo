classdef SliceViewer < handle
    % SliceViewer - Interactive slice visualization for 2D/3D data
    %
    % This viewer provides slice-by-slice visualization with interactive
    % controls for navigation, windowing, and analysis.
    %
    % Example:
    %   viewer = mufo.vis.SliceViewer('Parent', panel);
    %   viewer.setData(volume3D);
    
    properties (Access = private)
        parent              % Parent container
        panel               % Main panel
        ax                  % Axes handle
        imageHandle         % Image handle
        data                % Current data
        currentSlice        % Current slice index
        colormap            % Current colormap
        clim                % Color limits
        
        % UI Controls
        controls            % Control panel
        sliceSlider         % Slice navigation slider
        sliceEdit           % Slice number edit box
        windowSlider        % Window level slider
        levelSlider         % Window width slider
        colormapPopup       % Colormap selection
        
        % Interactive elements
        roiManager          % ROI management
        lineProfile         % Line profile tool
        pixelInfo           % Pixel information display
    end
    
    properties (Constant)
        DEFAULT_COLORMAP = 'gray';
        COLORMAPS = {'gray', 'bone', 'hot', 'jet', 'parula', 'viridis'};
    end
    
    events
        SliceChanged
        ROICreated
        WindowChanged
    end
    
    methods
        function obj = SliceViewer(varargin)
            % Constructor
            %
            % Parameters:
            %   'Parent' - Parent container (figure or panel)
            %   'Position' - Position vector [x y w h]
            
            p = inputParser;
            addParameter(p, 'Parent', [], @ishandle);
            addParameter(p, 'Position', [0 0 1 1], @isnumeric);
            parse(p, varargin{:});
            
            obj.parent = p.Results.Parent;
            if isempty(obj.parent)
                obj.parent = figure('Name', 'Slice Viewer');
            end
            
            % Create UI
            obj.createUI(p.Results.Position);
            
            % Initialize
            obj.colormap = obj.DEFAULT_COLORMAP;
            obj.currentSlice = 1;
        end
        
        function setData(obj, data)
            % Set data to visualize
            %
            % Input:
            %   data - 2D or 3D numeric array
            
            if ~isnumeric(data)
                error('SliceViewer:InvalidData', ...
                    'Data must be numeric array');
            end
            
            obj.data = data;
            
            % Reset view
            obj.currentSlice = 1;
            obj.clim = [min(data(:)), max(data(:))];
            
            % Update UI
            obj.updateSliceControls();
            obj.updateDisplay();
            obj.updateWindowControls();
            
            % Fit to view
            axis(obj.ax, 'image');
            
            % Notify
            notify(obj, 'SliceChanged');
        end
        
        function slice = getCurrentSlice(obj)
            % Get current slice data
            if ndims(obj.data) == 3
                slice = obj.data(:, :, obj.currentSlice);
            else
                slice = obj.data;
            end
        end
        
        function setSlice(obj, sliceNum)
            % Set current slice number
            if ndims(obj.data) == 3
                maxSlice = size(obj.data, 3);
                if sliceNum >= 1 && sliceNum <= maxSlice
                    obj.currentSlice = sliceNum;
                    obj.updateDisplay();
                    obj.updateSliceControls();
                    notify(obj, 'SliceChanged');
                end
            end
        end
        
        function setColormap(obj, cmapName)
            % Set colormap
            if any(strcmp(cmapName, obj.COLORMAPS))
                obj.colormap = cmapName;
                colormap(obj.ax, cmapName);
            end
        end
        
        function setWindow(obj, window, level)
            % Set window width and level
            %
            % Inputs:
            %   window - Window width
            %   level - Window center level
            
            obj.clim = [level - window/2, level + window/2];
            if ~isempty(obj.imageHandle) && isvalid(obj.imageHandle)
                obj.ax.CLim = obj.clim;
            end
            notify(obj, 'WindowChanged');
        end
        
        function roi = createROI(obj, type)
            % Create interactive ROI
            %
            % Input:
            %   type - ROI type ('rectangle', 'ellipse', 'polygon', 'line')
            
            switch lower(type)
                case 'rectangle'
                    roi = drawrectangle(obj.ax);
                case 'ellipse'
                    roi = drawellipse(obj.ax);
                case 'polygon'
                    roi = drawpolygon(obj.ax);
                case 'line'
                    roi = drawline(obj.ax);
                case 'freehand'
                    roi = drawfreehand(obj.ax);
                otherwise
                    error('SliceViewer:InvalidROI', ...
                        'Invalid ROI type: %s', type);
            end
            
            % Add listener for ROI changes
            addlistener(roi, 'MovingROI', @(~,~) obj.updateROIStats(roi));
            addlistener(roi, 'ROIMoved', @(~,~) obj.updateROIStats(roi));
            
            % Store ROI
            if isempty(obj.roiManager)
                obj.roiManager = roi;
            else
                obj.roiManager(end+1) = roi;
            end
            
            notify(obj, 'ROICreated');
            
            % Calculate initial stats
            obj.updateROIStats(roi);
        end
        
        function profile = showLineProfile(obj)
            % Show interactive line profile
            
            % Create line ROI
            lineROI = drawline(obj.ax);
            
            % Create profile figure
            profFig = figure('Name', 'Line Profile', ...
                'Position', [100, 100, 600, 400]);
            
            % Update profile function
            updateProfile = @(~,~) obj.plotProfile(profFig, lineROI);
            
            % Add listeners
            addlistener(lineROI, 'MovingROI', updateProfile);
            addlistener(lineROI, 'ROIMoved', updateProfile);
            addlistener(obj, 'SliceChanged', updateProfile);
            
            % Initial plot
            updateProfile();
            
            profile = struct('roi', lineROI, 'figure', profFig);
        end
        
        function show(obj)
            % Show the viewer
            obj.panel.Visible = 'on';
        end
        
        function hide(obj)
            % Hide the viewer
            obj.panel.Visible = 'off';
        end
        
        function exportImage(obj, filename)
            % Export current view as image
            if nargin < 2
                [file, path] = uiputfile({'*.png'; '*.tif'; '*.jpg'}, ...
                    'Export Image');
                if isequal(file, 0)
                    return;
                end
                filename = fullfile(path, file);
            end
            
            % Get current frame
            frame = getframe(obj.ax);
            imwrite(frame.cdata, filename);
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
            
            % Image axes (left side)
            obj.ax = axes(...
                'Parent', obj.panel, ...
                'Units', 'normalized', ...
                'Position', [0.02, 0.1, 0.65, 0.88], ...
                'Box', 'on', ...
                'XTick', [], ...
                'YTick', []);
            
            % Control panel (right side)
            obj.controls = uipanel(...
                'Parent', obj.panel, ...
                'Title', 'Controls', ...
                'Units', 'normalized', ...
                'Position', [0.68, 0.1, 0.30, 0.88]);
            
            % Create controls
            obj.createSliceControls();
            obj.createWindowControls();
            obj.createDisplayControls();
            obj.createROIControls();
            
            % Pixel info panel (bottom)
            obj.createPixelInfo();
            
            % Set up interactions
            obj.setupInteractions();
        end
        
        function createSliceControls(obj)
            % Create slice navigation controls
            
            yPos = 0.85;
            
            % Slice label
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Slice:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.2, 0.05], ...
                'HorizontalAlignment', 'left');
            
            % Slice slider
            obj.sliceSlider = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.06, 0.65, 0.05], ...
                'Min', 1, ...
                'Max', 1, ...
                'Value', 1, ...
                'SliderStep', [1, 10], ...
                'Callback', @(src,~) obj.onSliceSlider(src));
            
            % Slice edit box
            obj.sliceEdit = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'edit', ...
                'String', '1', ...
                'Units', 'normalized', ...
                'Position', [0.72, yPos-0.06, 0.23, 0.05], ...
                'Callback', @(src,~) obj.onSliceEdit(src));
            
            % Navigation buttons
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', '<<', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.12, 0.12, 0.05], ...
                'Callback', @(~,~) obj.firstSlice());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', '<', ...
                'Units', 'normalized', ...
                'Position', [0.19, yPos-0.12, 0.12, 0.05], ...
                'Callback', @(~,~) obj.previousSlice());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Play', ...
                'Units', 'normalized', ...
                'Position', [0.33, yPos-0.12, 0.29, 0.05], ...
                'Callback', @(src,~) obj.playSlices(src));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', '>', ...
                'Units', 'normalized', ...
                'Position', [0.64, yPos-0.12, 0.12, 0.05], ...
                'Callback', @(~,~) obj.nextSlice());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', '>>', ...
                'Units', 'normalized', ...
                'Position', [0.78, yPos-0.12, 0.12, 0.05], ...
                'Callback', @(~,~) obj.lastSlice());
        end
        
        function createWindowControls(obj)
            % Create window/level controls
            
            yPos = 0.65;
            
            % Window width
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Window:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.3, 0.05], ...
                'HorizontalAlignment', 'left');
            
            obj.windowSlider = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.06, 0.9, 0.05], ...
                'Min', 0, ...
                'Max', 1, ...
                'Value', 1, ...
                'Callback', @(~,~) obj.onWindowChange());
            
            % Window level
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Level:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.12, 0.3, 0.05], ...
                'HorizontalAlignment', 'left');
            
            obj.levelSlider = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.18, 0.9, 0.05], ...
                'Min', 0, ...
                'Max', 1, ...
                'Value', 0.5, ...
                'Callback', @(~,~) obj.onWindowChange());
            
            % Auto button
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Auto', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.24, 0.25, 0.05], ...
                'Callback', @(~,~) obj.autoWindow());
            
            % Reset button
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Reset', ...
                'Units', 'normalized', ...
                'Position', [0.35, yPos-0.24, 0.25, 0.05], ...
                'Callback', @(~,~) obj.resetWindow());
            
            % Presets button
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Presets', ...
                'Units', 'normalized', ...
                'Position', [0.65, yPos-0.24, 0.3, 0.05], ...
                'Callback', @(~,~) obj.showWindowPresets());
        end
        
        function createDisplayControls(obj)
            % Create display controls
            
            yPos = 0.35;
            
            % Colormap
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Colormap:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.3, 0.05], ...
                'HorizontalAlignment', 'left');
            
            obj.colormapPopup = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'popupmenu', ...
                'String', obj.COLORMAPS, ...
                'Units', 'normalized', ...
                'Position', [0.35, yPos, 0.6, 0.05], ...
                'Value', 1, ...
                'Callback', @(src,~) obj.onColormapChange(src));
            
            % Display options
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'checkbox', ...
                'String', 'Show Grid', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.08, 0.4, 0.05], ...
                'Callback', @(src,~) obj.toggleGrid(src));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'checkbox', ...
                'String', 'Show Colorbar', ...
                'Units', 'normalized', ...
                'Position', [0.5, yPos-0.08, 0.45, 0.05], ...
                'Callback', @(src,~) obj.toggleColorbar(src));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'checkbox', ...
                'String', 'Interpolate', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.14, 0.4, 0.05], ...
                'Value', 0, ...
                'Callback', @(src,~) obj.toggleInterpolation(src));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'checkbox', ...
                'String', 'Equal Aspect', ...
                'Units', 'normalized', ...
                'Position', [0.5, yPos-0.14, 0.45, 0.05], ...
                'Value', 1, ...
                'Callback', @(src,~) obj.toggleAspectRatio(src));
        end
        
        function createROIControls(obj)
            % Create ROI controls
            
            yPos = 0.15;
            
            % ROI tools label
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'ROI Tools:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.9, 0.05], ...
                'HorizontalAlignment', 'left', ...
                'FontWeight', 'bold');
            
            % ROI buttons
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Rectangle', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.06, 0.28, 0.05], ...
                'Callback', @(~,~) obj.createROI('rectangle'));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Ellipse', ...
                'Units', 'normalized', ...
                'Position', [0.36, yPos-0.06, 0.28, 0.05], ...
                'Callback', @(~,~) obj.createROI('ellipse'));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Polygon', ...
                'Units', 'normalized', ...
                'Position', [0.67, yPos-0.06, 0.28, 0.05], ...
                'Callback', @(~,~) obj.createROI('polygon'));
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Line Profile', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.12, 0.43, 0.05], ...
                'Callback', @(~,~) obj.showLineProfile());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Clear ROIs', ...
                'Units', 'normalized', ...
                'Position', [0.52, yPos-0.12, 0.43, 0.05], ...
                'Callback', @(~,~) obj.clearROIs());
        end
        
        function createPixelInfo(obj)
            % Create pixel information display
            obj.pixelInfo = uicontrol(...
                'Parent', obj.panel, ...
                'Style', 'text', ...
                'String', 'Position: (-, -) | Value: -', ...
                'Units', 'normalized', ...
                'Position', [0.02, 0.02, 0.65, 0.06], ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor', [0.95, 0.95, 0.95], ...
                'FontName', 'FixedWidth');
        end
        
        function setupInteractions(obj)
            % Set up mouse interactions
            set(obj.ax, 'ButtonDownFcn', @(~,~) obj.onAxesClick());
            
            % Mouse motion for pixel info
            set(obj.panel.Parent, 'WindowButtonMotionFcn', ...
                @(~,~) obj.updatePixelInfo());
        end
        
        function updateDisplay(obj)
            % Update image display
            if isempty(obj.data)
                return;
            end
            
            % Get current slice
            slice = obj.getCurrentSlice();
            
            % Display image
            if isempty(obj.imageHandle) || ~isvalid(obj.imageHandle)
                obj.imageHandle = imagesc(obj.ax, slice);
                axis(obj.ax, 'image');
                colormap(obj.ax, obj.colormap);
            else
                obj.imageHandle.CData = slice;
            end
            
            % Update color limits
            if ~isempty(obj.clim)
                obj.ax.CLim = obj.clim;
            end
            
            % Update title
            if ndims(obj.data) == 3
                title(obj.ax, sprintf('Slice %d / %d', ...
                    obj.currentSlice, size(obj.data, 3)));
            end
        end
        
        function updateSliceControls(obj)
            % Update slice navigation controls
            if ndims(obj.data) == 3
                maxSlice = size(obj.data, 3);
                
                % Update slider
                obj.sliceSlider.Max = maxSlice;
                obj.sliceSlider.Value = obj.currentSlice;
                
                if maxSlice > 1
                    obj.sliceSlider.SliderStep = [1/(maxSlice-1), 10/(maxSlice-1)];
                end
                
                % Update edit box
                obj.sliceEdit.String = num2str(obj.currentSlice);
                
                % Enable controls
                obj.sliceSlider.Enable = 'on';
                obj.sliceEdit.Enable = 'on';
            else
                % Disable for 2D data
                obj.sliceSlider.Enable = 'off';
                obj.sliceEdit.Enable = 'off';
            end
        end
        
        function updateWindowControls(obj)
            % Update window/level controls
            if isempty(obj.data)
                return;
            end
            
            dataMin = min(obj.data(:));
            dataMax = max(obj.data(:));
            dataRange = dataMax - dataMin;
            
            % Update window slider (width)
            obj.windowSlider.Min = 0;
            obj.windowSlider.Max = dataRange;
            obj.windowSlider.Value = obj.clim(2) - obj.clim(1);
            
            % Update level slider (center)
            obj.levelSlider.Min = dataMin;
            obj.levelSlider.Max = dataMax;
            obj.levelSlider.Value = mean(obj.clim);
        end
        
        function updatePixelInfo(obj)
            % Update pixel information display
            if isempty(obj.data) || isempty(obj.ax)
                return;
            end
            
            % Get current point
            cp = obj.ax.CurrentPoint;
            x = round(cp(1,1));
            y = round(cp(1,2));
            
            % Check bounds
            slice = obj.getCurrentSlice();
            [h, w] = size(slice);
            
            if x >= 1 && x <= w && y >= 1 && y <= h
                value = slice(y, x);
                
                % Format string
                if isinteger(value)
                    valStr = sprintf('%d', value);
                else
                    valStr = sprintf('%.4f', value);
                end
                
                infoStr = sprintf('Position: (%d, %d) | Value: %s', ...
                    x, y, valStr);
                
                % Add slice info for 3D
                if ndims(obj.data) == 3
                    infoStr = sprintf('%s | Slice: %d/%d', ...
                        infoStr, obj.currentSlice, size(obj.data, 3));
                end
                
                obj.pixelInfo.String = infoStr;
            end
        end
        
        function updateROIStats(obj, roi)
            % Update ROI statistics
            try
                slice = obj.getCurrentSlice();
                mask = roi.createMask(obj.imageHandle);
                
                % Calculate statistics
                values = slice(mask);
                
                stats = struct();
                stats.mean = mean(values);
                stats.std = std(values);
                stats.min = min(values);
                stats.max = max(values);
                stats.area = sum(mask(:));
                
                % Display stats
                statStr = sprintf('Mean: %.2f ± %.2f\nMin: %.2f, Max: %.2f\nArea: %d px', ...
                    stats.mean, stats.std, stats.min, stats.max, stats.area);
                
                % Update ROI label if it supports it
                if isprop(roi, 'Label')
                    roi.Label = statStr;
                end
            catch
                % Handle errors silently
            end
        end
        
        function plotProfile(obj, fig, lineROI)
            % Plot line profile in given figure
            figure(fig);
            clf;
            
            try
                % Get line endpoints
                pos = lineROI.Position;
                x = [pos(1,1), pos(2,1)];
                y = [pos(1,2), pos(2,2)];
                
                % Extract profile
                slice = obj.getCurrentSlice();
                profile = improfile(slice, x, y);
                
                % Plot
                subplot(2,1,1);
                plot(profile, 'b-', 'LineWidth', 2);
                xlabel('Distance (pixels)');
                ylabel('Intensity');
                title(sprintf('Line Profile - Slice %d', obj.currentSlice));
                grid on;
                
                % Show gradient
                subplot(2,1,2);
                grad = gradient(profile);
                plot(grad, 'r-', 'LineWidth', 2);
                xlabel('Distance (pixels)');
                ylabel('Gradient');
                title('Profile Gradient');
                grid on;
                
            catch
                % Handle errors
                clf;
                text(0.5, 0.5, 'Error plotting profile', ...
                    'HorizontalAlignment', 'center');
            end
        end
        
        % Callback functions
        function onSliceSlider(obj, src)
            % Handle slice slider change
            obj.currentSlice = round(src.Value);
            obj.updateDisplay();
            obj.sliceEdit.String = num2str(obj.currentSlice);
            notify(obj, 'SliceChanged');
        end
        
        function onSliceEdit(obj, src)
            % Handle slice edit box change
            val = str2double(src.String);
            if ~isnan(val) && val >= 1 && val <= size(obj.data, 3)
                obj.currentSlice = round(val);
                obj.updateDisplay();
                obj.sliceSlider.Value = obj.currentSlice;
                notify(obj, 'SliceChanged');
            else
                % Reset to current value
                src.String = num2str(obj.currentSlice);
            end
        end
        
        function onWindowChange(obj)
            % Handle window/level change
            window = obj.windowSlider.Value;
            level = obj.levelSlider.Value;
            obj.setWindow(window, level);
        end
        
        function onColormapChange(obj, src)
            % Handle colormap change
            cmapName = obj.COLORMAPS{src.Value};
            obj.setColormap(cmapName);
        end
        
        function onAxesClick(obj)
            % Handle axes click
            % Could be extended for ROI placement
        end
        
        % Navigation functions
        function firstSlice(obj)
            obj.setSlice(1);
        end
        
        function lastSlice(obj)
            if ndims(obj.data) == 3
                obj.setSlice(size(obj.data, 3));
            end
        end
        
        function previousSlice(obj)
            if obj.currentSlice > 1
                obj.setSlice(obj.currentSlice - 1);
            end
        end
        
        function nextSlice(obj)
            if ndims(obj.data) == 3 && obj.currentSlice < size(obj.data, 3)
                obj.setSlice(obj.currentSlice + 1);
            end
        end
        
        function playSlices(obj, btn)
            % Play through slices
            if strcmp(btn.String, 'Play')
                btn.String = 'Stop';
                
                % Animation loop
                while isvalid(btn) && strcmp(btn.String, 'Stop')
                    obj.nextSlice();
                    if obj.currentSlice == size(obj.data, 3)
                        obj.firstSlice();
                    end
                    pause(0.1);
                    drawnow;
                end
                
                if isvalid(btn)
                    btn.String = 'Play';
                end
            else
                btn.String = 'Play';
            end
        end
        
        % Window functions
        function autoWindow(obj)
            % Auto adjust window/level
            slice = obj.getCurrentSlice();
            
            % Use percentiles for robust auto-windowing
            vals = slice(:);
            vals(isnan(vals) | isinf(vals)) = [];
            
            p1 = prctile(vals, 1);
            p99 = prctile(vals, 99);
            
            obj.setWindow(p99 - p1, (p99 + p1) / 2);
            obj.updateWindowControls();
        end
        
        function resetWindow(obj)
            % Reset window to full data range
            if ~isempty(obj.data)
                obj.clim = [min(obj.data(:)), max(obj.data(:))];
                obj.updateDisplay();
                obj.updateWindowControls();
            end
        end
        
        function showWindowPresets(obj)
            % Show window presets menu
            presets = {
                'Soft Tissue', 400, 50;
                'Bone', 2000, 400;
                'Lung', 1500, -600;
                'Brain', 80, 40;
                'Liver', 150, 50;
                'Full Range', [], []
            };
            
            % Create menu
            menu = uicontextmenu(obj.panel.Parent);
            
            for i = 1:size(presets, 1)
                uimenu(menu, 'Text', presets{i,1}, ...
                    'Callback', @(~,~) obj.applyPreset(presets{i,2}, presets{i,3}));
            end
            
            % Show menu at current pointer location
            set(menu, 'Visible', 'on');
        end
        
        function applyPreset(obj, window, level)
            % Apply window preset
            if isempty(window)
                obj.resetWindow();
            else
                obj.setWindow(window, level);
                obj.updateWindowControls();
            end
        end
        
        % Display option functions
        function toggleGrid(obj, src)
            % Toggle grid display
            if src.Value
                grid(obj.ax, 'on');
            else
                grid(obj.ax, 'off');
            end
        end
        
        function toggleColorbar(obj, src)
            % Toggle colorbar
            if src.Value
                colorbar(obj.ax);
            else
                colorbar(obj.ax, 'off');
            end
        end
        
        function toggleInterpolation(obj, src)
            % Toggle image interpolation
            if ~isempty(obj.imageHandle) && isvalid(obj.imageHandle)
                if src.Value
                    obj.imageHandle.Interpolation = 'bilinear';
                else
                    obj.imageHandle.Interpolation = 'nearest';
                end
            end
        end
        
        function toggleAspectRatio(obj, src)
            % Toggle equal aspect ratio
            if src.Value
                axis(obj.ax, 'image');
            else
                axis(obj.ax, 'normal');
            end
        end
        
        function clearROIs(obj)
            % Clear all ROIs
            if ~isempty(obj.roiManager)
                for i = 1:length(obj.roiManager)
                    if isvalid(obj.roiManager(i))
                        delete(obj.roiManager(i));
                    end
                end
                obj.roiManager = [];
            end
        end
    end
end