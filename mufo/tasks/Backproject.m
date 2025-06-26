classdef Backproject < mufo.core.UFOTask
    % Backproject - UFO backprojection task wrapper
    %
    % Performs filtered backprojection reconstruction.
    % This is the final step in FBP reconstruction after filtering.
    %
    % Example:
    %   task = mufo.tasks.Backproject();
    %   task.setCenter(1024.5);              % Center of rotation
    %   task.setAngleRange(0, 2*pi, 1800);   % 0-360 degrees, 1800 projections
    
    properties (Access = private)
        backprojectInfo     % Information about backprojection
    end
    
    methods
        function obj = Backproject(varargin)
            % Constructor
            obj@mufo.core.UFOTask('backproject');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setCenter(obj, center)
            % Set center of rotation
            %
            % Input:
            %   center - Axis position (center of rotation in pixels)
            
            if center < 0
                error('Backproject:InvalidCenter', ...
                    'Center must be non-negative');
            end
            
            obj.setParameter('axis-pos', center);
        end
        
        function setAngleStep(obj, step)
            % Set angle step in radians
            %
            % Input:
            %   step - Angle increment between projections (radians)
            
            if step <= 0
                error('Backproject:InvalidAngleStep', ...
                    'Angle step must be positive');
            end
            
            obj.setParameter('angle-step', step);
        end
        
        function setAngleRange(obj, startAngle, endAngle, numProjections)
            % Set angle range and calculate step
            %
            % Inputs:
            %   startAngle - Starting angle in radians
            %   endAngle - Ending angle in radians
            %   numProjections - Number of projections (optional)
            
            obj.setParameter('angle-start', startAngle);
            obj.setParameter('angle-end', endAngle);
            
            % Calculate angle step if number of projections given
            if nargin > 3 && numProjections > 1
                angleRange = endAngle - startAngle;
                angleStep = angleRange / numProjections;
                obj.setAngleStep(angleStep);
            end
        end
        
        function setAngleRangeDegrees(obj, startDeg, endDeg, numProjections)
            % Set angle range in degrees
            %
            % Inputs:
            %   startDeg - Starting angle in degrees
            %   endDeg - Ending angle in degrees
            %   numProjections - Number of projections (optional)
            
            startRad = deg2rad(startDeg);
            endRad = deg2rad(endDeg);
            
            if nargin > 3
                obj.setAngleRange(startRad, endRad, numProjections);
            else
                obj.setAngleRange(startRad, endRad);
            end
        end
        
        function setPrecision(obj, mode)
            % Set precision mode
            %
            % Input:
            %   mode - 'single', 'half', or 'int8'
            
            validModes = {'single', 'half', 'int8'};
            if ~any(strcmpi(mode, validModes))
                error('Backproject:InvalidPrecision', ...
                    'Precision must be single, half, or int8');
            end
            
            obj.setParameter('precision-mode', mode);
        end
        
        function setAddressingMode(obj, mode)
            % Set addressing mode for out-of-bounds access
            %
            % Input:
            %   mode - 'clamp', 'clamp_to_edge', 'repeat', or 'mirrored_repeat'
            
            obj.setParameter('addressing-mode', mode);
        end
        
        function setOutputSize(obj, width, height)
            % Set output volume size
            %
            % Inputs:
            %   width - Output width
            %   height - Output height (optional, defaults to width)
            
            if nargin < 3
                height = width;  % Square output
            end
            
            obj.setParameter('output-width', width);
            obj.setParameter('output-height', height);
        end
        
        function info = calculateBackprojectionInfo(obj, sinogramDims)
            % Calculate backprojection parameters
            %
            % Input:
            %   sinogramDims - Sinogram dimensions [width, numAngles, height]
            % Output:
            %   info - Structure with backprojection information
            
            info = struct();
            info.sinogramDims = sinogramDims;
            
            % Get parameters
            info.center = obj.getParameter('axis-pos');
            info.angleStart = obj.getParameter('angle-start');
            info.angleEnd = obj.getParameter('angle-end');
            info.angleStep = obj.getParameter('angle-step');
            
            % Determine output size
            outputWidth = obj.getParameter('output-width');
            outputHeight = obj.getParameter('output-height');
            
            if isempty(outputWidth)
                % Default to detector width
                info.outputSize = [sinogramDims(1), sinogramDims(1)];
            else
                info.outputSize = [outputWidth, ...
                    ifelse(isempty(outputHeight), outputWidth, outputHeight)];
            end
            
            % Calculate angle information
            if ~isempty(info.angleStep)
                info.numAngles = round((info.angleEnd - info.angleStart) / info.angleStep);
            else
                info.numAngles = sinogramDims(2);
                if info.numAngles > 1
                    info.angleStep = (info.angleEnd - info.angleStart) / info.numAngles;
                end
            end
            
            % Memory estimate
            info.slicesPerChunk = sinogramDims(3);
            info.bytesPerVoxel = 4;  % Single precision
            info.memoryPerSlice = prod(info.outputSize) * info.bytesPerVoxel;
            info.totalMemory = info.memoryPerSlice * info.slicesPerChunk;
            
            obj.backprojectInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== Backprojection Parameters ===\n');
                fprintf('Sinogram: [%d × %d × %d]\n', sinogramDims);
                fprintf('Output volume: [%d × %d × %d]\n', ...
                    info.outputSize(1), info.outputSize(2), info.slicesPerChunk);
                fprintf('Center of rotation: %.2f\n', info.center);
                fprintf('Angle range: %.1f° to %.1f°\n', ...
                    rad2deg(info.angleStart), rad2deg(info.angleEnd));
                fprintf('Number of angles: %d\n', info.numAngles);
                fprintf('Angle step: %.4f° (%.6f rad)\n', ...
                    rad2deg(info.angleStep), info.angleStep);
                fprintf('Memory required: %.2f MB\n', info.totalMemory / 1e6);
            end
        end
        
        function center = findCenter(obj, sinogram, searchRange)
            % Find center of rotation using correlation method
            %
            % Inputs:
            %   sinogram - Single sinogram slice
            %   searchRange - Range to search around image center
            % Output:
            %   center - Estimated center of rotation
            
            if nargin < 3
                searchRange = 50;
            end
            
            [width, numAngles] = size(sinogram);
            imageCenter = width / 2;
            
            % Test range
            centers = imageCenter + (-searchRange:0.5:searchRange);
            correlation = zeros(size(centers));
            
            % Find opposite projections (0 and 180 degrees)
            idx1 = 1;
            idx2 = round(numAngles / 2);
            
            proj1 = sinogram(:, idx1);
            proj2 = sinogram(:, idx2);
            
            % Test each center
            for i = 1:length(centers)
                % Flip proj2 and shift according to center
                shift = 2 * (centers(i) - imageCenter);
                proj2_flipped = flipud(proj2);
                
                % Shift and compute correlation
                if shift > 0
                    proj2_shifted = [zeros(round(shift), 1); ...
                                    proj2_flipped(1:end-round(shift))];
                else
                    proj2_shifted = [proj2_flipped(-round(shift)+1:end); ...
                                    zeros(-round(shift), 1)];
                end
                
                % Compute correlation
                correlation(i) = corr(proj1, proj2_shifted);
            end
            
            % Find maximum correlation
            [~, maxIdx] = max(correlation);
            center = centers(maxIdx);
            
            % Plot if no output requested
            if nargout == 0
                figure('Name', 'Center of Rotation Finding');
                plot(centers, correlation, 'b-', 'LineWidth', 2);
                hold on;
                plot(center, correlation(maxIdx), 'ro', 'MarkerSize', 10);
                xlabel('Center Position (pixels)');
                ylabel('Correlation');
                title(sprintf('Center of Rotation: %.2f', center));
                grid on;
            end
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get reconstruction parameters
                reconParams = config.getReconstructionParams();
                
                % Set center if available
                if ~isempty(reconParams.center)
                    obj.setCenter(reconParams.center);
                end
                
                % Set angle range
                angleRange = reconParams.angleRange;
                if angleRange == 360
                    obj.setAngleRangeDegrees(0, 360);
                elseif angleRange == 180
                    obj.setAngleRangeDegrees(0, 180);
                else
                    obj.setAngleRangeDegrees(0, angleRange);
                end
                
                % Set precision
                if isfield(reconParams, 'precision')
                    obj.setPrecision(reconParams.precision);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'center')
                    obj.setCenter(config.center);
                end
                
                if isfield(config, 'angleStart') && isfield(config, 'angleEnd')
                    obj.setAngleRange(config.angleStart, config.angleEnd);
                end
                
                if isfield(config, 'angleStep')
                    obj.setAngleStep(config.angleStep);
                end
                
                if isfield(config, 'precision')
                    obj.setPrecision(config.precision);
                end
                
                if isfield(config, 'outputSize')
                    if isscalar(config.outputSize)
                        obj.setOutputSize(config.outputSize);
                    else
                        obj.setOutputSize(config.outputSize(1), config.outputSize(2));
                    end
                end
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define backprojection parameters
            
            % Center of rotation
            obj.addParameterDef('axis-pos', 'numeric', ...
                'description', 'Axis position (center of rotation)', ...
                'required', true);
            
            % Angle parameters
            obj.addParameterDef('angle-step', 'numeric', ...
                'description', 'Angle step in radians', ...
                'range', [0, inf]);
            
            obj.addParameterDef('angle-start', 'numeric', ...
                'description', 'Starting angle in radians', ...
                'default', 0);
            
            obj.addParameterDef('angle-end', 'numeric', ...
                'description', 'Ending angle in radians', ...
                'default', 2*pi);
            
            % Precision mode
            obj.addParameterDef('precision-mode', 'string', ...
                'description', 'Precision mode', ...
                'values', {'single', 'half', 'int8'}, ...
                'default', 'single');
            
            % Addressing mode
            obj.addParameterDef('addressing-mode', 'string', ...
                'description', 'Addressing mode for out-of-bounds', ...
                'values', {'clamp', 'clamp_to_edge', 'repeat', 'mirrored_repeat'}, ...
                'default', 'clamp_to_edge');
            
            % Output size
            obj.addParameterDef('output-width', 'integer', ...
                'description', 'Output width', ...
                'range', [1, inf]);
            
            obj.addParameterDef('output-height', 'integer', ...
                'description', 'Output height', ...
                'range', [1, inf]);
            
            % Volume parameters
            obj.addParameterDef('volume-width', 'integer', ...
                'description', 'Volume width', ...
                'range', [1, inf]);
            
            obj.addParameterDef('volume-height', 'integer', ...
                'description', 'Volume height', ...
                'range', [1, inf]);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            % Check if center is set
            if isempty(obj.getParameter('axis-pos'))
                valid = false;
                warning('Backproject:NoCenter', ...
                    'Center of rotation must be set');
            end
            
            % Check angle parameters
            angleStart = obj.getParameter('angle-start');
            angleEnd = obj.getParameter('angle-end');
            angleStep = obj.getParameter('angle-step');
            
            if angleEnd <= angleStart
                valid = false;
                warning('Backproject:InvalidAngleRange', ...
                    'End angle must be greater than start angle');
            end
            
            if ~isempty(angleStep) && angleStep <= 0
                valid = false;
                warning('Backproject:InvalidAngleStep', ...
                    'Angle step must be positive');
            end
            
            % Warn if angle step not set
            if isempty(angleStep)
                warning('Backproject:NoAngleStep', ...
                    'Angle step not set - UFO will use default');
            end
        end
    end
    
    methods (Access = private)
        function result = ifelse(obj, condition, trueVal, falseVal)
            % Helper function for conditional assignment
            if condition
                result = trueVal;
            else
                result = falseVal;
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate backprojection operations
            
            fprintf('\n=== Backprojection Task Demo ===\n\n');
            
            % Create task
            task = mufo.tasks.Backproject();
            
            % Example 1: Standard FBP backprojection
            fprintf('1. Standard FBP backprojection setup:\n');
            task.setCenter(1024.5);
            task.setAngleRangeDegrees(0, 360, 1800);
            task.calculateBackprojectionInfo([2048, 1800, 100]);
            
            % Example 2: Limited angle reconstruction
            fprintf('\n2. Limited angle (180°) reconstruction:\n');
            task.setAngleRangeDegrees(0, 180, 900);
            task.calculateBackprojectionInfo([2048, 900, 100]);
            
            % Example 3: High precision with custom output size
            fprintf('\n3. High precision with custom output:\n');
            task.setPrecision('single');
            task.setOutputSize(1024, 1024);
            task.calculateBackprojectionInfo([2048, 1800, 100]);
            
            % Example 4: Center finding demo
            fprintf('\n4. Center of rotation finding:\n');
            % Create synthetic sinogram
            angles = linspace(0, pi, 180);
            detector = linspace(-1, 1, 512);
            [D, A] = meshgrid(detector, angles);
            
            % Simple phantom (off-center circle)
            phantom = zeros(512, 180);
            for i = 1:180
                theta = angles(i);
                % Circle at (0.1, 0) with radius 0.3
                t = D(i,:) * cos(theta) + 0.1 * sin(theta);
                phantom(:,i) = abs(t) < 0.3;
            end
            
            % Add center offset
            trueCenter = 256 + 25;  % 25 pixels off
            phantom = circshift(phantom, [25, 0]);
            
            % Find center
            center = task.findCenter(phantom', 50);
            fprintf('True center: %.1f, Found center: %.1f\n', ...
                trueCenter, center);
            
            fprintf('\nDemo complete.\n');
        end
    end
end