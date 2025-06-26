classdef Crop < mufo.core.UFOTask
    % Crop - UFO crop task wrapper
    %
    % Crops data to a specified region of interest. Useful for removing
    % unwanted borders or focusing on specific areas.
    %
    % Example:
    %   task = mufo.tasks.Crop();
    %   task.setRegion(100, 100, 1848, 1848);  % x, y, width, height
    
    properties (Access = private)
        cropInfo        % Information about crop operation
    end
    
    methods
        function obj = Crop(varargin)
            % Constructor
            obj@mufo.core.UFOTask('crop');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setRegion(obj, x, y, width, height)
            % Set crop region
            %
            % Inputs:
            %   x - X offset (0-based)
            %   y - Y offset (0-based)
            %   width - Region width
            %   height - Region height
            
            if x < 0 || y < 0
                error('Crop:InvalidOffset', ...
                    'Offsets must be non-negative');
            end
            
            if width <= 0 || height <= 0
                error('Crop:InvalidSize', ...
                    'Width and height must be positive');
            end
            
            obj.setParameter('x', x);
            obj.setParameter('y', y);
            obj.setParameter('width', width);
            obj.setParameter('height', height);
        end
        
        function setCenterRegion(obj, inputWidth, inputHeight, cropWidth, cropHeight)
            % Set crop region centered in input
            %
            % Inputs:
            %   inputWidth - Input image width
            %   inputHeight - Input image height
            %   cropWidth - Desired crop width
            %   cropHeight - Desired crop height
            
            if cropWidth > inputWidth || cropHeight > inputHeight
                error('Crop:RegionTooLarge', ...
                    'Crop region exceeds input dimensions');
            end
            
            % Calculate centered offsets
            x = floor((inputWidth - cropWidth) / 2);
            y = floor((inputHeight - cropHeight) / 2);
            
            obj.setRegion(x, y, cropWidth, cropHeight);
        end
        
        function setSquareCrop(obj, inputWidth, inputHeight)
            % Set square crop from center of rectangular input
            %
            % Inputs:
            %   inputWidth - Input width
            %   inputHeight - Input height
            
            size = min(inputWidth, inputHeight);
            obj.setCenterRegion(inputWidth, inputHeight, size, size);
        end
        
        function info = getCropInfo(obj, inputDims)
            % Get information about crop operation
            %
            % Input:
            %   inputDims - Input dimensions [width, height] or [width, height, depth]
            % Output:
            %   info - Structure with crop information
            
            info = struct();
            info.inputDims = inputDims;
            
            % Get crop parameters
            info.x = obj.getParameter('x');
            info.y = obj.getParameter('y');
            info.width = obj.getParameter('width');
            info.height = obj.getParameter('height');
            
            % Calculate output dimensions
            outputDims = inputDims;
            outputDims(1) = info.width;
            outputDims(2) = info.height;
            info.outputDims = outputDims;
            
            % Check bounds
            info.xEnd = info.x + info.width;
            info.yEnd = info.y + info.height;
            
            if info.xEnd > inputDims(1) || info.yEnd > inputDims(2)
                warning('Crop:OutOfBounds', ...
                    'Crop region extends beyond input boundaries');
                info.valid = false;
            else
                info.valid = true;
            end
            
            % Calculate data reduction
            inputPixels = prod(inputDims(1:2));
            outputPixels = info.width * info.height;
            info.reductionFactor = inputPixels / outputPixels;
            info.percentRetained = 100 * outputPixels / inputPixels;
            
            obj.cropInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== Crop Operation ===\n');
                fprintf('Input: [%s]\n', num2str(inputDims));
                fprintf('Crop region: (%d, %d) to (%d, %d)\n', ...
                    info.x, info.y, info.xEnd, info.yEnd);
                fprintf('Output: [%s]\n', num2str(info.outputDims));
                fprintf('Data retained: %.1f%%\n', info.percentRetained);
                fprintf('Reduction factor: %.2fx\n', info.reductionFactor);
                if ~info.valid
                    fprintf('WARNING: Crop region out of bounds!\n');
                end
            end
        end
        
        function region = detectContentRegion(obj, data, threshold)
            % Detect region containing actual content
            %
            % Inputs:
            %   data - 2D image data
            %   threshold - Threshold for content detection
            % Output:
            %   region - Struct with x, y, width, height
            
            if nargin < 3
                threshold = 0.01 * max(data(:));
            end
            
            % Find non-background pixels
            mask = data > threshold;
            
            % Find bounding box
            [rows, cols] = find(mask);
            if isempty(rows)
                warning('Crop:NoContent', 'No content found above threshold');
                region = struct('x', 0, 'y', 0, ...
                               'width', size(data, 2), ...
                               'height', size(data, 1));
                return;
            end
            
            minRow = min(rows);
            maxRow = max(rows);
            minCol = min(cols);
            maxCol = max(cols);
            
            % Add small margin
            margin = 10;
            minRow = max(1, minRow - margin);
            maxRow = min(size(data, 1), maxRow + margin);
            minCol = max(1, minCol - margin);
            maxCol = min(size(data, 2), maxCol + margin);
            
            % Convert to crop parameters (0-based for UFO)
            region = struct();
            region.x = minCol - 1;
            region.y = minRow - 1;
            region.width = maxCol - minCol + 1;
            region.height = maxRow - minRow + 1;
            
            % Display result
            figure('Name', 'Content Detection');
            
            subplot(1,2,1);
            imagesc(data);
            colormap gray;
            hold on;
            rectangle('Position', [minCol, minRow, region.width, region.height], ...
                     'EdgeColor', 'r', 'LineWidth', 2);
            title('Detected Content Region');
            
            subplot(1,2,2);
            imagesc(data(minRow:maxRow, minCol:maxCol));
            colormap gray;
            title('Cropped Result');
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Check for crop settings
                cropConfig = config.get('preprocessing.crop', struct());
                
                if isfield(cropConfig, 'x') && isfield(cropConfig, 'y') && ...
                   isfield(cropConfig, 'width') && isfield(cropConfig, 'height')
                    obj.setRegion(cropConfig.x, cropConfig.y, ...
                                 cropConfig.width, cropConfig.height);
                elseif isfield(cropConfig, 'mode')
                    % Handle special crop modes
                    if strcmpi(cropConfig.mode, 'square')
                        % Will need input dimensions at runtime
                        warning('Crop:NeedInputDims', ...
                            'Square crop mode requires input dimensions');
                    end
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'x') && isfield(config, 'y') && ...
                   isfield(config, 'width') && isfield(config, 'height')
                    obj.setRegion(config.x, config.y, ...
                                 config.width, config.height);
                end
                
                if isfield(config, 'center') && config.center
                    % Centered crop - needs input dimensions
                    if isfield(config, 'inputWidth') && isfield(config, 'inputHeight')
                        obj.setCenterRegion(config.inputWidth, config.inputHeight, ...
                                          config.width, config.height);
                    end
                end
            end
        end
        
        function visualizeCrop(obj, data)
            % Visualize crop operation on sample data
            %
            % Input:
            %   data - 2D image to demonstrate cropping
            
            [height, width] = size(data);
            
            % Get crop region
            x = obj.getParameter('x');
            y = obj.getParameter('y');
            cropWidth = obj.getParameter('width');
            cropHeight = obj.getParameter('height');
            
            % Perform crop (1-based for MATLAB)
            cropped = data(y+1:y+cropHeight, x+1:x+cropWidth);
            
            % Create visualization
            figure('Name', 'Crop Visualization');
            
            % Original with overlay
            subplot(2,2,1);
            imagesc(data);
            colormap gray;
            hold on;
            rectangle('Position', [x+1, y+1, cropWidth, cropHeight], ...
                     'EdgeColor', 'r', 'LineWidth', 2);
            title('Original with Crop Region');
            xlabel(sprintf('Size: %d × %d', width, height));
            
            % Cropped result
            subplot(2,2,2);
            imagesc(cropped);
            colormap gray;
            title('Cropped Result');
            xlabel(sprintf('Size: %d × %d', cropWidth, cropHeight));
            
            % Profile comparison
            subplot(2,2,3);
            midY = round(height/2);
            midYCrop = round(cropHeight/2);
            plot(data(midY, :), 'b-', 'LineWidth', 1.5);
            hold on;
            plot(x+1:x+cropWidth, cropped(midYCrop, :), 'r-', 'LineWidth', 2);
            xlabel('X Position');
            ylabel('Intensity');
            title('Horizontal Profiles');
            legend('Original', 'Cropped', 'Location', 'best');
            grid on;
            
            % Statistics
            subplot(2,2,4);
            axis off;
            text(0.1, 0.9, sprintf('Original size: %d × %d', width, height));
            text(0.1, 0.8, sprintf('Cropped size: %d × %d', cropWidth, cropHeight));
            text(0.1, 0.7, sprintf('Offset: (%d, %d)', x, y));
            text(0.1, 0.6, sprintf('Data retained: %.1f%%', ...
                100 * cropWidth * cropHeight / (width * height)));
            
            text(0.1, 0.4, 'Statistics:', 'FontWeight', 'bold');
            text(0.1, 0.3, sprintf('Original mean: %.3f', mean(data(:))));
            text(0.1, 0.2, sprintf('Cropped mean: %.3f', mean(cropped(:))));
            text(0.1, 0.1, sprintf('Original std: %.3f', std(data(:))));
            text(0.1, 0.0, sprintf('Cropped std: %.3f', std(cropped(:))));
            
            title('Crop Information');
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define crop parameters
            
            % Position
            obj.addParameterDef('x', 'integer', ...
                'description', 'X offset (0-based)', ...
                'range', [0, inf], ...
                'required', true);
            task
            obj.addParameterDef('y', 'integer', ...
                'description', 'Y offset (0-based)', ...
                'range', [0, inf], ...
                'required', true);
            
            % Size
            obj.addParameterDef('width', 'integer', ...
                'description', 'Crop width', ...
                'range', [1, inf], ...
                'required', true);
            
            obj.addParameterDef('height', 'integer', ...
                'description', 'Crop height', ...
                'range', [1, inf], ...
                'required', true);
            
            % Options
            obj.addParameterDef('out-of-bounds-mode', 'string', ...
                'description', 'Behavior for out-of-bounds regions', ...
                'values', {'error', 'clamp', 'pad'}, ...
                'default', 'error');
            
            obj.addParameterDef('pad-value', 'numeric', ...
                'description', 'Padding value for out-of-bounds', ...
                'default', 0);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            % Check all required parameters are set
            x = obj.getParameter('x');
            y = obj.getParameter('y');
            width = obj.getParameter('width');
            height = obj.getParameter('height');
            
            if isempty(x) || isempty(y) || isempty(width) || isempty(height)
                valid = false;
                warning('Crop:MissingParameters', ...
                    'All crop parameters (x, y, width, height) must be set');
                return;
            end
            task
            % Check for valid values
            if x < 0 || y < 0
                valid = false;
                warning('Crop:NegativeOffset', ...
                    'Crop offsets must be non-negative');
            end
            
            if width <= 0 || height <= 0
                valid = false;
                warning('Crop:InvalidSize', ...
                    'Crop width and height must be positive');
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate crop operations
            
            fprintf('\n=== Crop Task Demo ===\n\n');
            
            % Create synthetic image with border artifacts
            baseImage = phantom(512);
            
            % Add border artifacts
            image = baseImage;
            image(1:50, :) = 0.2;      % Top border
            image(end-49:end, :) = 0.2; % Bottom border
            image(:, 1:50) = 0.2;       % Left border
            image(:, end-49:end) = 0.2; % Right border
            
            % Add some noise
            image = image + 0.05 * randn(size(image));
            
            %task Create task
            task = mufo.tasks.Crop();
            
            % Example 1: Manual crop to remove borders
            fprintf('1. Manual crop to remove borders:\n');
            task.setRegion(50, 50, 412, 412);
            task.getCropInfo([512, 512]);
            
            % Example 2: Center crop to square
            fprintf('\n2. Center crop to square:\n');
            task.setSquareCrop(640, 480);
            task.getCropInfo([640, 480]);
            
            % Example 3: Automatic content detection
            fprintf('\n3. Automatic content detection:\n');
            region = task.detectContentRegion(image, 0.3);
            task.setRegion(region.x, region.y, region.width, region.height);
            
            % Example 4: Visualize crop
            fprintf('\n4. Visualizing crop operation:\n');
            task.visualizeCrop(image);
            
            fprintf('\nDemo complete.\n');
        end
    end
end