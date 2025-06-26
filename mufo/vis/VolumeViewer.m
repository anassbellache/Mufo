classdef VolumeViewer < handle
    % VolumeViewer - 3D volume visualization with multiple rendering modes
    %
    % This viewer provides various 3D visualization methods including
    % orthogonal slices, MIP, and volume rendering.
    %
    % Example:
    %   viewer = mufo.vis.VolumeViewer('Parent', panel);
    %   viewer.setData(volume3D);
    %   viewer.setRenderMode('mip');
    
    properties (Access = private)
        parent              % Parent container
        panel               % Main panel
        ax                  % Main axes (for 3D view)
        sliceAxes           % Axes for orthogonal slices
        data                % Volume data
        
        % Rendering
        renderMode          % Current render mode
        volumeHandle        % Volume graphics handle
        sliceHandles        % Slice graphics handles
        isoHandle           % Isosurface handle
        
        % UI Controls
        controls            % Control panel
        modeSelector        % Render mode selector
        opacitySlider       % Opacity control
        thresholdSlider     % Threshold control
        colorMapSelector    % Colormap selector
        
        % Slice positions
        sliceX              % X slice position
        sliceY              % Y slice position
        sliceZ              % Z slice position
        
        % Display parameters
        colormap            % Current colormap
        alphamap            % Current alphamap
        clim                % Color limits
        threshold           % Isosurface threshold
    end
    
    properties (Constant)
        RENDER_MODES = {'Orthogonal Slices', 'MIP', 'Volume Render', ...
                       'Isosurface', 'Multi-Planar'};
        COLORMAPS = {'gray', 'bone', 'hot', 'jet', 'parula', 'viridis'};
        DEFAULT_COLORMAP = 'gray';
    end
    
    events
        ViewChanged
        SliceChanged
    end
    
    methods
        function obj = VolumeViewer(varargin)
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
                obj.parent = figure('Name', 'Volume Viewer');
            end
            
            % Create UI
            obj.createUI(p.Results.Position);
            
            % Initialize
            obj.renderMode = 'Orthogonal Slices';
            obj.colormap = obj.DEFAULT_COLORMAP;
        end
        
        function setData(obj, data, varargin)
            % Set volume data
            %
            % Inputs:
            %   data - 3D volume data
            %   'Spacing' - Voxel spacing [x y z]
            
            p = inputParser;
            addRequired(p, 'data', @(x) isnumeric(x) && ndims(x) == 3);
            addParameter(p, 'Spacing', [1 1 1], @isnumeric);
            parse(p, data, varargin{:});
            
            obj.data = data;
            
            % Initialize slice positions at center
            [ny, nx, nz] = size(data);
            obj.sliceX = round(nx / 2);
            obj.sliceY = round(ny / 2);
            obj.sliceZ = round(nz / 2);
            
            % Set color limits
            obj.clim = [min(data(:)), max(data(:))];
            
            % Set default threshold
            obj.threshold = mean(obj.clim);
            
            % Update display
            obj.updateDisplay();
            obj.updateControls();
            
            % Notify
            notify(obj, 'ViewChanged');
        end
        
        function setRenderMode(obj, mode)
            % Set rendering mode
            %
            % Input:
            %   mode - Render mode name
            
            if any(strcmp(mode, obj.RENDER_MODES))
                obj.renderMode = mode;
                obj.updateDisplay();
                notify(obj, 'ViewChanged');
            else
                error('VolumeViewer:InvalidMode', ...
                    'Invalid render mode: %s', mode);
            end
        end
        
        function setSlicePosition(obj, dim, pos)
            % Set slice position
            %
            % Inputs:
            %   dim - Dimension ('x', 'y', or 'z')
            %   pos - Slice position
            
            switch lower(dim)
                case 'x'
                    if pos >= 1 && pos <= size(obj.data, 2)
                        obj.sliceX = round(pos);
                    end
                case 'y'
                    if pos >= 1 && pos <= size(obj.data, 1)
                        obj.sliceY = round(pos);
                    end
                case 'z'
                    if pos >= 1 && pos <= size(obj.data, 3)
                        obj.sliceZ = round(pos);
                    end
            end
            
            obj.updateSlices();
            notify(obj, 'SliceChanged');
        end
        
        function createIsosurface(obj, value)
            % Create isosurface at specific value
            %
            % Input:
            %   value - Isosurface value
            
            if nargin < 2
                value = obj.threshold;
            end
            
            % Clear existing isosurface
            if ~isempty(obj.isoHandle) && isvalid(obj.isoHandle)
                delete(obj.isoHandle);
            end
            
            % Create isosurface
            obj.isoHandle = patch(isosurface(obj.data, value), ...
                'Parent', obj.ax, ...
                'FaceColor', [0.5, 0.5, 0.8], ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.7);
            
            % Add lighting
            camlight(obj.ax);
            lighting(obj.ax, 'gouraud');
            
            % Update view
            view(obj.ax, 3);
            axis(obj.ax, 'vis3d');
            rotate3d(obj.ax, 'on');
        end
        
        function mip = createMIP(obj, direction)
            % Create Maximum Intensity Projection
            %
            % Input:
            %   direction - Projection direction ('x', 'y', 'z', or 'all')
            % Output:
            %   mip - MIP image(s)
            
            if nargin < 2
                direction = 'all';
            end
            
            switch lower(direction)
                case 'x'
                    mip = squeeze(max(obj.data, [], 2));
                case 'y'
                    mip = squeeze(max(obj.data, [], 1));
                case 'z'
                    mip = squeeze(max(obj.data, [], 3));
                case 'all'
                    mip = struct();
                    mip.x = squeeze(max(obj.data, [], 2));
                    mip.y = squeeze(max(obj.data, [], 1));
                    mip.z = squeeze(max(obj.data, [], 3));
            end
        end
        
        function show(obj)
            % Show the viewer
            obj.panel.Visible = 'on';
        end
        
        function hide(obj)
            % Hide the viewer
            obj.panel.Visible = 'off';
        end
        
        function exportView(obj, filename)
            % Export current view
            if nargin < 2
                [file, path] = uiputfile({'*.png'; '*.jpg'; '*.tif'}, ...
                    'Export View');
                if isequal(file, 0)
                    return;
                end
                filename = fullfile(path, file);
            end
            
            % Get current view
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
            
            % Create layout based on initial mode
            obj.createOrthogonalLayout();
            
            % Control panel (right side)
            obj.controls = uipanel(...
                'Parent', obj.panel, ...
                'Title', 'Volume Controls', ...
                'Units', 'normalized', ...
                'Position', [0.75, 0.02, 0.24, 0.96]);
            
            % Create controls
            obj.createVolumeControls();
        end
        
        function createOrthogonalLayout(obj)
            % Create layout for orthogonal slices
            
            % Clear any existing axes
            obj.clearAxes();
            
            % Main 3D view (top left)
            obj.ax = axes(...
                'Parent', obj.panel, ...
                'Units', 'normalized', ...
                'Position', [0.02, 0.51, 0.35, 0.45], ...
                'Box', 'on');
            title(obj.ax, '3D View');
            
            % XY slice (top right)
            obj.sliceAxes.xy = axes(...
                'Parent', obj.panel, ...
                'Units', 'normalized', ...
                'Position', [0.39, 0.51, 0.35, 0.45], ...
                'Box', 'on');
            title(obj.sliceAxes.xy, 'XY Slice (Axial)');
            
            % XZ slice (bottom left)
            obj.sliceAxes.xz = axes(...
                'Parent', obj.panel, ...
                'Units', 'normalized', ...
                'Position', [0.02, 0.02, 0.35, 0.45], ...
                'Box', 'on');
            title(obj.sliceAxes.xz, 'XZ Slice (Coronal)');
            
            % YZ slice (bottom right)
            obj.sliceAxes.yz = axes(...
                'Parent', obj.panel, ...
                'Units', 'normalized', ...
                'Position', [0.39, 0.02, 0.35, 0.45], ...
                'Box', 'on');
            title(obj.sliceAxes.yz, 'YZ Slice (Sagittal)');
            
            % Set up interactions
            obj.setupSliceInteractions();
        end
        
        function createSingleLayout(obj)
            % Create layout for single 3D view
            
            % Clear any existing axes
            obj.clearAxes();
            
            % Single large 3D view
            obj.ax = axes(...
                'Parent', obj.panel, ...
                'Units', 'normalized', ...
                'Position', [0.02, 0.02, 0.72, 0.94], ...
                'Box', 'on');
            
            % Enable 3D rotation
            rotate3d(obj.ax, 'on');
        end
        
        function createVolumeControls(obj)
            % Create volume control UI
            
            yPos = 0.9;
            
            % Render mode
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Render Mode:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.9, 0.05], ...
                'HorizontalAlignment', 'left');
            
            obj.modeSelector = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'popupmenu', ...
                'String', obj.RENDER_MODES, ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.06, 0.9, 0.05], ...
                'Callback', @(src,~) obj.onModeChange(src));
            
            % Colormap
            yPos = yPos - 0.15;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Colormap:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.9, 0.05], ...
                'HorizontalAlignment', 'left');
            
            obj.colorMapSelector = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'popupmenu', ...
                'String', obj.COLORMAPS, ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.06, 0.9, 0.05], ...
                'Callback', @(src,~) obj.onColormapChange(src));
            
            % Opacity control
            yPos = yPos - 0.15;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Opacity:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.9, 0.05], ...
                'HorizontalAlignment', 'left');
            
            obj.opacitySlider = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.06, 0.9, 0.05], ...
                'Min', 0, 'Max', 1, 'Value', 0.5, ...
                'Callback', @(~,~) obj.updateOpacity());
            
            % Threshold control
            yPos = yPos - 0.15;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Threshold:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.9, 0.05], ...
                'HorizontalAlignment', 'left');
            
            obj.thresholdSlider = uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos-0.06, 0.9, 0.05], ...
                'Min', 0, 'Max', 1, 'Value', 0.5, ...
                'Callback', @(~,~) obj.updateThreshold());
            
            % Slice controls
            yPos = yPos - 0.15;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Slice Controls:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.9, 0.05], ...
                'HorizontalAlignment', 'left', ...
                'FontWeight', 'bold');
            
            % X slice
            yPos = yPos - 0.08;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'X:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.1, 0.05], ...
                'HorizontalAlignment', 'left');
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.15, yPos, 0.7, 0.05], ...
                'Tag', 'SliceX', ...
                'Callback', @(src,~) obj.onSliceChange(src, 'x'));
            
            % Y slice
            yPos = yPos - 0.08;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Y:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.1, 0.05], ...
                'HorizontalAlignment', 'left');
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.15, yPos, 0.7, 0.05], ...
                'Tag', 'SliceY', ...
                'Callback', @(src,~) obj.onSliceChange(src, 'y'));
            
            % Z slice
            yPos = yPos - 0.08;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'text', ...
                'String', 'Z:', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.1, 0.05], ...
                'HorizontalAlignment', 'left');
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'slider', ...
                'Units', 'normalized', ...
                'Position', [0.15, yPos, 0.7, 0.05], ...
                'Tag', 'SliceZ', ...
                'Callback', @(src,~) obj.onSliceChange(src, 'z'));
            
            % Action buttons
            yPos = yPos - 0.12;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Reset View', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.42, 0.08], ...
                'Callback', @(~,~) obj.resetView());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Auto Contrast', ...
                'Units', 'normalized', ...
                'Position', [0.53, yPos, 0.42, 0.08], ...
                'Callback', @(~,~) obj.autoContrast());
            
            % Export options
            yPos = yPos - 0.1;
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Export View', ...
                'Units', 'normalized', ...
                'Position', [0.05, yPos, 0.42, 0.08], ...
                'Callback', @(~,~) obj.exportView());
            
            uicontrol(...
                'Parent', obj.controls, ...
                'Style', 'pushbutton', ...
                'String', 'Export Video', ...
                'Units', 'normalized', ...
                'Position', [0.53, yPos, 0.42, 0.08], ...
                'Callback', @(~,~) obj.exportRotation());
            
            % Volume info
            yPos = 0.15;
            infoPanel = uipanel(...
                'Parent', obj.controls, ...
                'Title', 'Volume Info', ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.02, 0.9, yPos-0.02]);
            
            obj.createInfoDisplay(infoPanel);
        end
        
        function createInfoDisplay(obj, parent)
            % Create volume information display
            
            if isempty(obj.data)
                return;
            end
            
            [ny, nx, nz] = size(obj.data);
            
            info = sprintf(['Size: %d x %d x %d\n', ...
                           'Voxels: %.1f M\n', ...
                           'Memory: %.1f MB\n', ...
                           'Range: [%.2f, %.2f]'], ...
                           nx, ny, nz, ...
                           numel(obj.data) / 1e6, ...
                           numel(obj.data) * 4 / 1e6, ...
                           min(obj.data(:)), max(obj.data(:)));
            
            uicontrol(...
                'Parent', parent, ...
                'Style', 'text', ...
                'String', info, ...
                'Units', 'normalized', ...
                'Position', [0.05, 0.05, 0.9, 0.9], ...
                'HorizontalAlignment', 'left', ...
                'FontName', 'FixedWidth');
        end
        
        function updateDisplay(obj)
            % Update display based on render mode
            
            if isempty(obj.data)
                return;
            end
            
            % Clear existing graphics
            obj.clearGraphics();
            
            switch obj.renderMode
                case 'Orthogonal Slices'
                    obj.createOrthogonalLayout();
                    obj.displayOrthogonalSlices();
                    
                case 'MIP'
                    obj.createSingleLayout();
                    obj.displayMIP();
                    
                case 'Volume Render'
                    obj.createSingleLayout();
                    obj.displayVolumeRender();
                    
                case 'Isosurface'
                    obj.createSingleLayout();
                    obj.displayIsosurface();
                    
                case 'Multi-Planar'
                    obj.createSingleLayout();
                    obj.displayMultiPlanar();
            end
            
            % Update colormap
            colormap(obj.ax, obj.colormap);
        end
        
        function displayOrthogonalSlices(obj)
            % Display orthogonal slice views
            
            % 3D view with slice planes
            axes(obj.ax);
            hold on;
            
            % Create slice planes
            [ny, nx, nz] = size(obj.data);
            
            % XY plane
            [X, Y] = meshgrid(1:nx, 1:ny);
            Z = ones(size(X)) * obj.sliceZ;
            obj.sliceHandles.xy = surf(X, Y, Z, squeeze(obj.data(:, :, obj.sliceZ)), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.8);
            
            % XZ plane
            [X, Z] = meshgrid(1:nx, 1:nz);
            Y = ones(size(X)) * obj.sliceY;
            obj.sliceHandles.xz = surf(X, Y, Z, squeeze(obj.data(obj.sliceY, :, :))', ...
                'EdgeColor', 'none', 'FaceAlpha', 0.8);
            
            % YZ plane
            [Y, Z] = meshgrid(1:ny, 1:nz);
            X = ones(size(Y)) * obj.sliceX;
            obj.sliceHandles.yz = surf(X, Y, Z, squeeze(obj.data(:, obj.sliceX, :))', ...
                'EdgeColor', 'none', 'FaceAlpha', 0.8);
            
            % Set properties
            axis(obj.ax, 'vis3d');
            xlabel(obj.ax, 'X');
            ylabel(obj.ax, 'Y');
            zlabel(obj.ax, 'Z');
            view(obj.ax, 3);
            rotate3d(obj.ax, 'on');
            
            % Update individual slice views
            obj.updateSlices();
        end
        
        function updateSlices(obj)
            % Update slice displays
            
            if isempty(obj.data) || ~strcmp(obj.renderMode, 'Orthogonal Slices')
                return;
            end
            
            % XY slice
            if isfield(obj.sliceAxes, 'xy') && isvalid(obj.sliceAxes.xy)
                axes(obj.sliceAxes.xy);
                imagesc(squeeze(obj.data(:, :, obj.sliceZ)));
                axis image;
                colormap(obj.colormap);
                caxis(obj.clim);
                title(sprintf('XY Slice (Z = %d)', obj.sliceZ));
                
                % Add crosshairs
                hold on;
                plot([obj.sliceX, obj.sliceX], ylim, 'r-', 'LineWidth', 1);
                plot(xlim, [obj.sliceY, obj.sliceY], 'g-', 'LineWidth', 1);
                hold off;
            end
            
            % XZ slice
            if isfield(obj.sliceAxes, 'xz') && isvalid(obj.sliceAxes.xz)
                axes(obj.sliceAxes.xz);
                imagesc(squeeze(obj.data(obj.sliceY, :, :))');
                axis image;
                colormap(obj.colormap);
                caxis(obj.clim);
                title(sprintf('XZ Slice (Y = %d)', obj.sliceY));
                
                % Add crosshairs
                hold on;
                plot([obj.sliceX, obj.sliceX], ylim, 'r-', 'LineWidth', 1);
                plot(xlim, [obj.sliceZ, obj.sliceZ], 'b-', 'LineWidth', 1);
                hold off;
            end
            
            % YZ slice
            if isfield(obj.sliceAxes, 'yz') && isvalid(obj.sliceAxes.yz)
                axes(obj.sliceAxes.yz);
                imagesc(squeeze(obj.data(:, obj.sliceX, :))');
                axis image;
                colormap(obj.colormap);
                caxis(obj.clim);
                title(sprintf('YZ Slice (X = %d)', obj.sliceX));
                
                % Add crosshairs
                hold on;
                plot([obj.sliceY, obj.sliceY], ylim, 'g-', 'LineWidth', 1);
                plot(xlim, [obj.sliceZ, obj.sliceZ], 'b-', 'LineWidth', 1);
                hold off;
            end
            
            % Update 3D slice positions if they exist
            if strcmp(obj.renderMode, 'Orthogonal Slices')
                obj.update3DSlices();
            end
        end
        
        function update3DSlices(obj)
            % Update 3D slice positions
            
            if isempty(obj.sliceHandles)
                return;
            end
            
            [ny, nx, nz] = size(obj.data);
            
            % Update XY slice
            if isfield(obj.sliceHandles, 'xy') && isvalid(obj.sliceHandles.xy)
                obj.sliceHandles.xy.ZData = ones(ny, nx) * obj.sliceZ;
                obj.sliceHandles.xy.CData = squeeze(obj.data(:, :, obj.sliceZ));
            end
            
            % Update XZ slice
            if isfield(obj.sliceHandles, 'xz') && isvalid(obj.sliceHandles.xz)
                obj.sliceHandles.xz.YData = ones(nz, nx) * obj.sliceY;
                obj.sliceHandles.xz.CData = squeeze(obj.data(obj.sliceY, :, :))';
            end
            
            % Update YZ slice
            if isfield(obj.sliceHandles, 'yz') && isvalid(obj.sliceHandles.yz)
                obj.sliceHandles.yz.XData = ones(ny, nz) * obj.sliceX;
                obj.sliceHandles.yz.CData = squeeze(obj.data(:, obj.sliceX, :))';
            end
        end
        
        function displayMIP(obj)
            % Display Maximum Intensity Projection
            
            % Create MIPs
            mips = obj.createMIP('all');
            
            % Display in subplots
            subplot(2,2,1, 'Parent', obj.ax.Parent);
            imagesc(mips.z);
            axis image;
            colormap(obj.colormap);
            caxis(obj.clim);
            title('MIP Z (Top View)');
            colorbar;
            
            subplot(2,2,2, 'Parent', obj.ax.Parent);
            imagesc(mips.y');
            axis image;
            colormap(obj.colormap);
            caxis(obj.clim);
            title('MIP Y (Front View)');
            colorbar;
            
            subplot(2,2,3, 'Parent', obj.ax.Parent);
            imagesc(mips.x');
            axis image;
            colormap(obj.colormap);
            caxis(obj.clim);
            title('MIP X (Side View)');
            colorbar;
            
            % 3D MIP visualization
            subplot(2,2,4, 'Parent', obj.ax.Parent);
            obj.ax = gca;
            
            % Create 3D box with MIPs on faces
            [ny, nx, nz] = size(obj.data);
            
            % Box vertices
            vertices = [1 1 1; nx 1 1; nx ny 1; 1 ny 1; ...
                       1 1 nz; nx 1 nz; nx ny nz; 1 ny nz];
            
            % Box edges
            edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; ...
                    1 5; 2 6; 3 7; 4 8];
            
            % Draw box
            hold on;
            for i = 1:size(edges, 1)
                plot3([vertices(edges(i,1),1), vertices(edges(i,2),1)], ...
                      [vertices(edges(i,1),2), vertices(edges(i,2),2)], ...
                      [vertices(edges(i,1),3), vertices(edges(i,2),3)], ...
                      'k-', 'LineWidth', 2);
            end
            
            view(3);
            axis vis3d;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('3D MIP Box');
            rotate3d on;
        end
        
        function displayVolumeRender(obj)
            % Display volume rendering
            
            % Note: MATLAB's volume rendering capabilities are limited
            % This is a simplified implementation
            
            % Create volumetric object
            % For actual volume rendering, consider using:
            % - volshow (requires Image Processing Toolbox)
            % - custom ray casting implementation
            % - External tools like VTK
            
            % Simplified approach using isocaps
            axes(obj.ax);
            
            % Multiple isosurfaces for depth
            nLevels = 5;
            levels = linspace(obj.clim(1), obj.clim(2), nLevels+2);
            levels = levels(2:end-1);
            
            alphaValues = linspace(0.1, 0.3, nLevels);
            colors = jet(nLevels);
            
            hold on;
            for i = 1:nLevels
                try
                    p = patch(isosurface(obj.data, levels(i)));
                    p.FaceColor = colors(i,:);
                    p.EdgeColor = 'none';
                    p.FaceAlpha = alphaValues(i);
                catch
                    % Handle case where isosurface fails
                end
            end
            
            % Add lighting
            camlight;
            lighting gouraud;
            
            % Set view
            view(3);
            axis vis3d;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Volume Rendering');
            rotate3d on;
            
            % Add note about limitations
            text(0.02, 0.98, 'Note: Simplified volume rendering', ...
                'Units', 'normalized', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', 'w', ...
                'EdgeColor', 'k');
        end
        
        function displayIsosurface(obj)
            % Display isosurface rendering
            
            axes(obj.ax);
            
            % Create isosurface
            obj.createIsosurface(obj.threshold);
            
            % Add bounding box
            [ny, nx, nz] = size(obj.data);
            
            hold on;
            % Draw box edges
            plot3([1 nx], [1 1], [1 1], 'k-');
            plot3([1 nx], [ny ny], [1 1], 'k-');
            plot3([1 nx], [1 1], [nz nz], 'k-');
            plot3([1 nx], [ny ny], [nz nz], 'k-');
            
            plot3([1 1], [1 ny], [1 1], 'k-');
            plot3([nx nx], [1 ny], [1 1], 'k-');
            plot3([1 1], [1 ny], [nz nz], 'k-');
            plot3([nx nx], [1 ny], [nz nz], 'k-');
            
            plot3([1 1], [1 1], [1 nz], 'k-');
            plot3([nx nx], [1 1], [1 nz], 'k-');
            plot3([1 1], [ny ny], [1 nz], 'k-');
            plot3([nx nx], [ny ny], [1 nz], 'k-');
            hold off;
            
            % Labels
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title(sprintf('Isosurface (Threshold = %.2f)', obj.threshold));
        end
        
        function displayMultiPlanar(obj)
            % Display multiple planar reconstructions
            
            % Create custom MPR view
            axes(obj.ax);
            
            % Show multiple slices along each axis
            nSlices = 5;
            
            [ny, nx, nz] = size(obj.data);
            
            hold on;
            
            % Z slices
            zSlices = round(linspace(1, nz, nSlices));
            for i = 1:length(zSlices)
                z = zSlices(i);
                [X, Y] = meshgrid(1:nx, 1:ny);
                Z = ones(size(X)) * z;
                
                h = surf(X, Y, Z, squeeze(obj.data(:, :, z)));
                h.EdgeColor = 'none';
                h.FaceAlpha = 0.3;
            end
            
            % Set properties
            axis vis3d;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            view(3);
            rotate3d on;
            colormap(obj.colormap);
            caxis(obj.clim);
            title('Multi-Planar Reconstruction');
        end
        
        function setupSliceInteractions(obj)
            % Set up mouse interactions for slice selection
            
            if isfield(obj.sliceAxes, 'xy')
                set(obj.sliceAxes.xy, 'ButtonDownFcn', ...
                    @(ax,evt) obj.onSliceClick(ax, evt, 'xy'));
            end
            
            if isfield(obj.sliceAxes, 'xz')
                set(obj.sliceAxes.xz, 'ButtonDownFcn', ...
                    @(ax,evt) obj.onSliceClick(ax, evt, 'xz'));
            end
            
            if isfield(obj.sliceAxes, 'yz')
                set(obj.sliceAxes.yz, 'ButtonDownFcn', ...
                    @(ax,evt) obj.onSliceClick(ax, evt, 'yz'));
            end
        end
        
        function updateControls(obj)
            % Update control values based on data
            
            if isempty(obj.data)
                return;
            end
            
            [ny, nx, nz] = size(obj.data);
            
            % Update slice sliders
            xSlider = findobj(obj.controls, 'Tag', 'SliceX');
            if ~isempty(xSlider)
                xSlider.Min = 1;
                xSlider.Max = nx;
                xSlider.Value = obj.sliceX;
                xSlider.SliderStep = [1/(nx-1), 10/(nx-1)];
            end
            
            ySlider = findobj(obj.controls, 'Tag', 'SliceY');
            if ~isempty(ySlider)
                ySlider.Min = 1;
                ySlider.Max = ny;
                ySlider.Value = obj.sliceY;
                ySlider.SliderStep = [1/(ny-1), 10/(ny-1)];
            end
            
            zSlider = findobj(obj.controls, 'Tag', 'SliceZ');
            if ~isempty(zSlider)
                zSlider.Min = 1;
                zSlider.Max = nz;
                zSlider.Value = obj.sliceZ;
                zSlider.SliderStep = [1/(nz-1), 10/(nz-1)];
            end
            
            % Update threshold slider
            if ~isempty(obj.thresholdSlider)
                obj.thresholdSlider.Min = obj.clim(1);
                obj.thresholdSlider.Max = obj.clim(2);
                obj.thresholdSlider.Value = obj.threshold;
            end
        end
        
        function clearAxes(obj)
            % Clear all axes
            if ~isempty(obj.ax) && isvalid(obj.ax)
                delete(obj.ax);
            end
            
            if ~isempty(obj.sliceAxes)
                fields = fieldnames(obj.sliceAxes);
                for i = 1:length(fields)
                    if isvalid(obj.sliceAxes.(fields{i}))
                        delete(obj.sliceAxes.(fields{i}));
                    end
                end
            end
            
            obj.sliceAxes = struct();
        end
        
        function clearGraphics(obj)
            % Clear graphics objects
            if ~isempty(obj.volumeHandle) && isvalid(obj.volumeHandle)
                delete(obj.volumeHandle);
            end
            
            if ~isempty(obj.isoHandle) && isvalid(obj.isoHandle)
                delete(obj.isoHandle);
            end
            
            if ~isempty(obj.sliceHandles)
                fields = fieldnames(obj.sliceHandles);
                for i = 1:length(fields)
                    if isvalid(obj.sliceHandles.(fields{i}))
                        delete(obj.sliceHandles.(fields{i}));
                    end
                end
            end
            
            obj.volumeHandle = [];
            obj.isoHandle = [];
            obj.sliceHandles = struct();
        end
        
        function resetView(obj)
            % Reset view to default
            if ~isempty(obj.ax) && isvalid(obj.ax)
                view(obj.ax, 3);
                axis(obj.ax, 'vis3d');
            end
            
            % Reset slice positions
            if ~isempty(obj.data)
                [ny, nx, nz] = size(obj.data);
                obj.sliceX = round(nx / 2);
                obj.sliceY = round(ny / 2);
                obj.sliceZ = round(nz / 2);
                
                obj.updateSlices();
                obj.updateControls();
            end
        end
        
        function autoContrast(obj)
            % Auto adjust contrast
            if isempty(obj.data)
                return;
            end
            
            % Use percentiles for robust contrast
            vals = obj.data(:);
            vals(isnan(vals) | isinf(vals)) = [];
            
            p1 = prctile(vals, 1);
            p99 = prctile(vals, 99);
            
            obj.clim = [p1, p99];
            
            % Update displays
            obj.updateDisplay();
        end
        
        function exportRotation(obj)
            % Export rotating view as video
            
            [file, path] = uiputfile({'*.mp4'; '*.avi'}, 'Export Video');
            if isequal(file, 0)
                return;
            end
            
            filename = fullfile(path, file);
            
            % Create video writer
            v = VideoWriter(filename);
            v.FrameRate = 30;
            open(v);
            
            % Rotate and capture frames
            nFrames = 360;
            
            for i = 1:nFrames
                view(obj.ax, [i, 30]);
                drawnow;
                
                frame = getframe(obj.ax);
                writeVideo(v, frame);
            end
            
            close(v);
            
            msgbox(sprintf('Video exported to: %s', filename), 'Export Complete');
        end
        
        % Callback functions
        function onModeChange(obj, src)
            % Handle render mode change
            modes = obj.RENDER_MODES;
            obj.setRenderMode(modes{src.Value});
        end
        
        function onColormapChange(obj, src)
            % Handle colormap change
            cmaps = obj.COLORMAPS;
            obj.colormap = cmaps{src.Value};
            
            % Apply to all axes
            if ~isempty(obj.ax) && isvalid(obj.ax)
                colormap(obj.ax, obj.colormap);
            end
            
            if ~isempty(obj.sliceAxes)
                fields = fieldnames(obj.sliceAxes);
                for i = 1:length(fields)
                    if isvalid(obj.sliceAxes.(fields{i}))
                        colormap(obj.sliceAxes.(fields{i}), obj.colormap);
                    end
                end
            end
        end
        
        function updateOpacity(obj)
            % Update opacity/transparency
            alpha = obj.opacitySlider.Value;
            
            % Apply to current graphics
            if strcmp(obj.renderMode, 'Isosurface') && ~isempty(obj.isoHandle)
                obj.isoHandle.FaceAlpha = alpha;
            elseif strcmp(obj.renderMode, 'Orthogonal Slices') && ~isempty(obj.sliceHandles)
                fields = fieldnames(obj.sliceHandles);
                for i = 1:length(fields)
                    if isvalid(obj.sliceHandles.(fields{i}))
                        obj.sliceHandles.(fields{i}).FaceAlpha = alpha;
                    end
                end
            end
        end
        
        function updateThreshold(obj)
            % Update threshold value
            range = obj.clim(2) - obj.clim(1);
            obj.threshold = obj.clim(1) + obj.thresholdSlider.Value * range;
            
            % Update isosurface if active
            if strcmp(obj.renderMode, 'Isosurface')
                obj.createIsosurface(obj.threshold);
            end
        end
        
        function onSliceChange(obj, src, dim)
            % Handle slice position change
            pos = round(src.Value);
            obj.setSlicePosition(dim, pos);
        end
        
        function onSliceClick(obj, ax, evt, plane)
            % Handle click on slice view
            pt = evt.IntersectionPoint;
            
            switch plane
                case 'xy'
                    obj.sliceX = round(pt(1));
                    obj.sliceY = round(pt(2));
                case 'xz'
                    obj.sliceX = round(pt(1));
                    obj.sliceZ = round(pt(2));
                case 'yz'
                    obj.sliceY = round(pt(1));
                    obj.sliceZ = round(pt(2));
            end
            
            obj.updateSlices();
            obj.updateControls();
            notify(obj, 'SliceChanged');
        end
    end
end