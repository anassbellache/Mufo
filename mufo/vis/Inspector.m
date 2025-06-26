classdef Inspector < handle
    % Inspector - Main inspection GUI for MUFO pipeline visualization
    %
    % This class provides the main interface for inspecting data at various
    % stages of the tomographic reconstruction pipeline.
    %
    % Example:
    %   inspector = mufo.vis.Inspector();
    %   inspector.inspect(data, 'Flat-field Corrected');
    
    properties (Access = private)
        fig                 % Main figure handle
        panels              % UI panels
        currentData         % Current dataset being inspected
        currentStage        % Current pipeline stage name
        viewers             % Specialized viewer instances
        toolbar             % Custom toolbar
        statusBar           % Status bar handle
        dataInfo            % Information about current data
        history             % Inspection history
    end
    
    properties (Constant)
        WINDOW_SIZE = [100, 100, 1600, 900];
        COLORS = struct(...
            'bg', [0.94, 0.94, 0.94], ...
            'panel', [1, 1, 1], ...
            'accent', [0.2, 0.4, 0.8]);
    end
    
    events
        DataChanged
        ViewChanged
        ExportRequested
    end
    
    methods
        function obj = Inspector(varargin)
            % Constructor - Initialize the inspection GUI
            %
            % Optional parameters:
            %   'Theme' - Color theme ('light', 'dark')
            %   'AutoDetect' - Auto-detect data type (default: true)
            
            p = inputParser;
            addParameter(p, 'Theme', 'light', @ischar);
            addParameter(p, 'AutoDetect', true, @islogical);
            parse(p, varargin{:});
            
            % Initialize components
            obj.panels = struct();
            obj.viewers = struct();
            obj.history = {};
            
            % Create GUI
            obj.createMainWindow();
            obj.createPanels();
            obj.createToolbar();
            obj.createMenus();
            obj.createStatusBar();
            
            % Initialize viewers
            obj.initializeViewers();
        end
        
        function inspect(obj, data, stageName, varargin)
            % Inspect data from a specific pipeline stage
            %
            % Inputs:
            %   data - Data to inspect (2D/3D array)
            %   stageName - Name of the pipeline stage
            %   'Metadata' - Additional metadata struct
            
            p = inputParser;
            addRequired(p, 'data', @isnumeric);
            addRequired(p, 'stageName', @ischar);
            addParameter(p, 'Metadata', struct(), @isstruct);
            parse(p, data, stageName, varargin{:});
            
            % Store data
            obj.currentData = data;
            obj.currentStage = stageName;
            
            % Analyze data
            obj.dataInfo = obj.analyzeData(data);
            obj.dataInfo.metadata = p.Results.Metadata;
            
            % Update UI
            obj.updateTitle(stageName);
            obj.updateStatusBar(sprintf('Loaded: %s', stageName));
            
            % Auto-detect and display appropriate view
            obj.autoDisplayData();
            
            % Add to history
            obj.addToHistory(stageName, obj.dataInfo);
            
            % Notify listeners
            notify(obj, 'DataChanged');
        end
        
        function show(obj)
            % Show the inspector window
            if ~isempty(obj.fig) && isvalid(obj.fig)
                figure(obj.fig);
            else
                obj.createMainWindow();
            end
        end
        
        function hide(obj)
            % Hide the inspector window
            if ~isempty(obj.fig) && isvalid(obj.fig)
                obj.fig.Visible = 'off';
            end
        end
        
        function exportData(obj, format, filename)
            % Export current view or data
            %
            % Inputs:
            %   format - Export format ('png', 'tiff', 'mat', 'h5')
            %   filename - Output filename
            
            if nargin < 3
                [file, path] = uiputfile(...
                    {'*.png', 'PNG Image'; ...
                     '*.tif', 'TIFF Stack'; ...
                     '*.mat', 'MATLAB Data'; ...
                     '*.h5', 'HDF5 File'}, ...
                    'Export Data');
                if isequal(file, 0)
                    return;
                end
                filename = fullfile(path, file);
                [~, ~, ext] = fileparts(filename);
                format = ext(2:end);
            end
            
            switch lower(format)
                case 'png'
                    obj.exportImage(filename);
                case {'tif', 'tiff'}
                    obj.exportTiff(filename);
                case 'mat'
                    obj.exportMat(filename);
                case {'h5', 'hdf5'}
                    obj.exportHDF5(filename);
                otherwise
                    error('Inspector:UnsupportedFormat', ...
                        'Unsupported export format: %s', format);
            end
            
            obj.updateStatusBar(sprintf('Exported to: %s', filename));
        end
    end
    
    methods (Access = private)
        function createMainWindow(obj)
            % Create the main inspector window
            obj.fig = figure(...
                'Name', 'MUFO Data Inspector', ...
                'NumberTitle', 'off', ...
                'Position', obj.WINDOW_SIZE, ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'Color', obj.COLORS.bg, ...
                'CloseRequestFcn', @(~,~) obj.onClose());
            
            % Set resize behavior
            obj.fig.SizeChangedFcn = @(~,~) obj.onResize();
        end
        
        function createPanels(obj)
            % Create UI panels
            
            % Main viewing panel (left side)
            obj.panels.main = uipanel(...
                'Parent', obj.fig, ...
                'Title', 'Data View', ...
                'Position', [0.01, 0.15, 0.68, 0.83], ...
                'BackgroundColor', obj.COLORS.panel, ...
                'BorderType', 'line', ...
                'HighlightColor', obj.COLORS.accent);
            
            % Control panel (right side)
            obj.panels.control = uipanel(...
                'Parent', obj.fig, ...
                'Title', 'Controls', ...
                'Position', [0.70, 0.15, 0.29, 0.83], ...
                'BackgroundColor', obj.COLORS.panel);
            
            % Info panel (bottom)
            obj.panels.info = uipanel(...
                'Parent', obj.fig, ...
                'Title', 'Data Information', ...
                'Position', [0.01, 0.01, 0.98, 0.13], ...
                'BackgroundColor', obj.COLORS.panel);
        end
        
        function createToolbar(obj)
            % Create custom toolbar
            obj.toolbar = uitoolbar(obj.fig);
            
            % Load icon data (use built-in icons or create custom ones)
            iconPath = fullfile(matlabroot, 'toolbox', 'matlab', 'icons');
            
            % File operations
            uipushtool(obj.toolbar, ...
                'CData', obj.loadIcon('file_open.png'), ...
                'TooltipString', 'Open Data', ...
                'ClickedCallback', @(~,~) obj.openData());
            
            uipushtool(obj.toolbar, ...
                'CData', obj.loadIcon('file_save.png'), ...
                'TooltipString', 'Export Data', ...
                'ClickedCallback', @(~,~) obj.exportData());
            
            % Add separator
            uipushtool(obj.toolbar, 'Separator', 'on', 'Visible', 'off');
            
            % View controls
            uitoggletool(obj.toolbar, ...
                'CData', obj.loadIcon('tool_zoom_in.png'), ...
                'TooltipString', 'Zoom', ...
                'ClickedCallback', @(src,~) obj.toggleZoom(src));
            
            uitoggletool(obj.toolbar, ...
                'CData', obj.loadIcon('tool_pan.png'), ...
                'TooltipString', 'Pan', ...
                'ClickedCallback', @(src,~) obj.togglePan(src));
            
            uitoggletool(obj.toolbar, ...
                'CData', obj.loadIcon('tool_datacursor.png'), ...
                'TooltipString', 'Data Cursor', ...
                'ClickedCallback', @(src,~) obj.toggleDataCursor(src));
            
            % Analysis tools
            uipushtool(obj.toolbar, ...
                'CData', obj.loadIcon('tool_line.png'), ...
                'TooltipString', 'Line Profile', ...
                'Separator', 'on', ...
                'ClickedCallback', @(~,~) obj.showLineProfile());
            
            uipushtool(obj.toolbar, ...
                'CData', obj.loadIcon('histogram.png'), ...
                'TooltipString', 'Histogram', ...
                'ClickedCallback', @(~,~) obj.showHistogram());
        end
        
        function createMenus(obj)
            % Create menu bar
            
            % File menu
            fileMenu = uimenu(obj.fig, 'Text', '&File');
            uimenu(fileMenu, 'Text', '&Open...', ...
                'Accelerator', 'O', ...
                'Callback', @(~,~) obj.openData());
            uimenu(fileMenu, 'Text', '&Export...', ...
                'Accelerator', 'S', ...
                'Callback', @(~,~) obj.exportData());
            uimenu(fileMenu, 'Text', 'Export &Batch...', ...
                'Separator', 'on', ...
                'Callback', @(~,~) obj.exportBatch());
            uimenu(fileMenu, 'Text', 'E&xit', ...
                'Separator', 'on', ...
                'Accelerator', 'Q', ...
                'Callback', @(~,~) obj.onClose());
            
            % View menu
            viewMenu = uimenu(obj.fig, 'Text', '&View');
            uimenu(viewMenu, 'Text', '&Slice View', ...
                'Callback', @(~,~) obj.switchToSliceView());
            uimenu(viewMenu, 'Text', '&Sinogram View', ...
                'Callback', @(~,~) obj.switchToSinogramView());
            uimenu(viewMenu, 'Text', '&3D Volume', ...
                'Callback', @(~,~) obj.switchToVolumeView());
            uimenu(viewMenu, 'Text', '&Comparison', ...
                'Separator', 'on', ...
                'Callback', @(~,~) obj.openComparison());
            
            % Tools menu
            toolsMenu = uimenu(obj.fig, 'Text', '&Tools');
            uimenu(toolsMenu, 'Text', '&Histogram', ...
                'Callback', @(~,~) obj.showHistogram());
            uimenu(toolsMenu, 'Text', '&Line Profile', ...
                'Callback', @(~,~) obj.showLineProfile());
            uimenu(toolsMenu, 'Text', '&ROI Analysis', ...
                'Callback', @(~,~) obj.showROIAnalysis());
            uimenu(toolsMenu, 'Text', '&Metrics', ...
                'Separator', 'on', ...
                'Callback', @(~,~) obj.showMetrics());
            
            % Window menu
            windowMenu = uimenu(obj.fig, 'Text', '&Window');
            uimenu(windowMenu, 'Text', '&Reset Layout', ...
                'Callback', @(~,~) obj.resetLayout());
            uimenu(windowMenu, 'Text', '&Full Screen', ...
                'Accelerator', 'F', ...
                'Callback', @(~,~) obj.toggleFullScreen());
        end
        
        function createStatusBar(obj)
            % Create status bar at bottom
            obj.statusBar = uicontrol(...
                'Parent', obj.fig, ...
                'Style', 'text', ...
                'String', 'Ready', ...
                'HorizontalAlignment', 'left', ...
                'Units', 'normalized', ...
                'Position', [0, 0, 1, 0.025], ...
                'BackgroundColor', [0.9, 0.9, 0.9]);
        end
        
        function initializeViewers(obj)
            % Initialize specialized viewer instances
            obj.viewers.slice = mufo.vis.SliceViewer('Parent', obj.panels.main);
            obj.viewers.sinogram = mufo.vis.SinogramViewer('Parent', obj.panels.main);
            obj.viewers.volume = mufo.vis.VolumeViewer('Parent', obj.panels.main);
            obj.viewers.metrics = mufo.vis.MetricsDisplay('Parent', obj.panels.info);
            
            % Hide all initially
            obj.viewers.slice.hide();
            obj.viewers.sinogram.hide();
            obj.viewers.volume.hide();
        end
        
        function info = analyzeData(obj, data)
            % Analyze data characteristics
            info = struct();
            info.size = size(data);
            info.ndims = ndims(data);
            info.class = class(data);
            info.min = min(data(:));
            info.max = max(data(:));
            info.mean = mean(data(:));
            info.std = std(data(:));
            info.memory = numel(data) * obj.sizeof(class(data)) / 1e6; % MB
            
            % Detect data type
            if info.ndims == 2
                info.type = 'slice';
            elseif info.ndims == 3
                % Check if it's a sinogram or volume
                if obj.isSinogram(data)
                    info.type = 'sinogram';
                else
                    info.type = 'volume';
                end
            else
                info.type = 'unknown';
            end
            
            % Check for special values
            info.hasNaN = any(isnan(data(:)));
            info.hasInf = any(isinf(data(:)));
            info.hasNegative = any(data(:) < 0);
        end
        
        function autoDisplayData(obj)
            % Automatically display data based on type
            
            % Hide all viewers first
            obj.hideAllViewers();
            
            % Show appropriate viewer
            switch obj.dataInfo.type
                case 'slice'
                    obj.viewers.slice.show();
                    obj.viewers.slice.setData(obj.currentData);
                    
                case 'sinogram'
                    obj.viewers.sinogram.show();
                    obj.viewers.sinogram.setData(obj.currentData);
                    
                case 'volume'
                    if obj.dataInfo.memory > 1000 % > 1GB
                        answer = questdlg(...
                            sprintf('Large volume detected (%.1f GB). Continue with 3D view?', ...
                                obj.dataInfo.memory/1000), ...
                            'Large Volume', 'Yes', 'Show Slices', 'Cancel', 'Show Slices');
                        
                        if strcmp(answer, 'Show Slices')
                            obj.viewers.slice.show();
                            obj.viewers.slice.setData(obj.currentData);
                            return;
                        elseif strcmp(answer, 'Cancel')
                            return;
                        end
                    end
                    
                    obj.viewers.volume.show();
                    obj.viewers.volume.setData(obj.currentData);
                    
                otherwise
                    warning('Inspector:UnknownDataType', ...
                        'Unknown data type. Showing as slices.');
                    obj.viewers.slice.show();
                    obj.viewers.slice.setData(obj.currentData);
            end
            
            % Update metrics
            obj.viewers.metrics.updateMetrics(obj.dataInfo);
        end
        
        function tf = isSinogram(obj, data)
            % Heuristic to detect if data is a sinogram
            tf = false;
            
            if ndims(data) == 3
                % Check aspect ratio and characteristics
                [ny, nx, nz] = size(data);
                
                % Sinograms typically have:
                % - Second dimension (detector) ≈ third dimension (slices)
                % - First dimension (angles) is different
                if abs(nx - nz) < 0.1 * nx && ny ~= nx
                    tf = true;
                end
                
                % Additional check: look for sinusoidal patterns
                midSlice = squeeze(data(:, round(nx/2), round(nz/2)));
                fftData = abs(fft(midSlice));
                
                % Sinograms have characteristic frequency content
                % This is a simplified check
                if max(fftData(2:10)) > 2 * mean(fftData(11:end))
                    tf = true;
                end
            end
        end
        
        function hideAllViewers(obj)
            % Hide all viewer panels
            viewers = fieldnames(obj.viewers);
            for i = 1:length(viewers)
                if ~strcmp(viewers{i}, 'metrics') % Keep metrics visible
                    obj.viewers.(viewers{i}).hide();
                end
            end
        end
        
        function updateTitle(obj, stageName)
            % Update window title
            if ~isempty(obj.dataInfo)
                sizeStr = sprintf('%dx%dx%d', obj.dataInfo.size);
                obj.fig.Name = sprintf('MUFO Inspector - %s [%s]', stageName, sizeStr);
            else
                obj.fig.Name = sprintf('MUFO Inspector - %s', stageName);
            end
        end
        
        function updateStatusBar(obj, message)
            % Update status bar message
            obj.statusBar.String = sprintf(' %s | %s', ...
                datestr(now, 'HH:MM:SS'), message);
        end
        
        function addToHistory(obj, stageName, info)
            % Add inspection to history
            entry = struct();
            entry.stage = stageName;
            entry.timestamp = now;
            entry.info = info;
            
            obj.history{end+1} = entry;
            
            % Limit history size
            if length(obj.history) > 50
                obj.history(1) = [];
            end
        end
        
        function icon = loadIcon(obj, filename)
            % Load icon for toolbar
            try
                iconPath = fullfile(matlabroot, 'toolbox', 'matlab', 'icons', filename);
                if exist(iconPath, 'file')
                    [icon, ~, alpha] = imread(iconPath);
                    if ~isempty(alpha)
                        % Handle transparency
                        icon = double(icon) / 255;
                        alpha = double(alpha) / 255;
                        bgColor = obj.COLORS.bg;
                        for i = 1:3
                            icon(:,:,i) = icon(:,:,i) .* alpha + ...
                                         bgColor(i) .* (1 - alpha);
                        end
                    end
                    % Resize to 16x16 if needed
                    if size(icon, 1) ~= 16 || size(icon, 2) ~= 16
                        icon = imresize(icon, [16, 16]);
                    end
                else
                    % Create default icon
                    icon = ones(16, 16, 3) * 0.5;
                end
            catch
                % Fallback to default icon
                icon = ones(16, 16, 3) * 0.5;
            end
        end
        
        function bytes = sizeof(obj, datatype)
            % Get size in bytes for data type
            switch datatype
                case {'single', 'int32', 'uint32'}
                    bytes = 4;
                case {'double', 'int64', 'uint64'}
                    bytes = 8;
                case {'int16', 'uint16'}
                    bytes = 2;
                case {'int8', 'uint8', 'logical'}
                    bytes = 1;
                otherwise
                    bytes = 8; % Default to double
            end
        end
        
        function onClose(obj)
            % Handle window close
            answer = questdlg('Close MUFO Inspector?', ...
                'Close Window', 'Yes', 'No', 'No');
            
            if strcmp(answer, 'Yes')
                delete(obj.fig);
            end
        end
        
        function onResize(obj)
            % Handle window resize
            % Update panel positions if needed
        end
        
        % Tool callbacks
        function toggleZoom(obj, src)
            % Toggle zoom mode
            if src.Value
                zoom(obj.fig, 'on');
                % Turn off other modes
                pan(obj.fig, 'off');
                datacursormode(obj.fig, 'off');
            else
                zoom(obj.fig, 'off');
            end
        end
        
        function togglePan(obj, src)
            % Toggle pan mode
            if src.Value
                pan(obj.fig, 'on');
                % Turn off other modes
                zoom(obj.fig, 'off');
                datacursormode(obj.fig, 'off');
            else
                pan(obj.fig, 'off');
            end
        end
        
        function toggleDataCursor(obj, src)
            % Toggle data cursor mode
            if src.Value
                datacursormode(obj.fig, 'on');
                % Turn off other modes
                zoom(obj.fig, 'off');
                pan(obj.fig, 'off');
            else
                datacursormode(obj.fig, 'off');
            end
        end
    end
end