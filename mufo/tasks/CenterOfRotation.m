classdef CenterOfRotation < mufo.core.UFOTask
    % CenterOfRotation - UFO center of rotation detection task wrapper
    %
    % Uses UFO's native correlate-stacks filter to find the center of
    % rotation by correlating projections at 0 and 180 degrees.
    %
    % The algorithm correlates opposite projections to find the shift
    % that maximizes correlation, which corresponds to the center of rotation.
    %
    % Example:
    %   task = mufo.tasks.CenterOfRotation();
    %   task.setSearchRange(50);        % Search ±50 pixels around center
    %   task.setSearchStep(0.5);        % Sub-pixel accuracy
    %   
    %   % Create workflow
    %   chain = mufo.core.UFOChain();
    %   chain.addTask(readTask);
    %   chain.addTask(task);
    %   chain.execute();
    %
    % Note: This task requires projections covering at least 180 degrees
    
    properties (Access = private)
        correlationInfo     % Information about correlation operation
        centerResult        % Detected center of rotation
    end
    
    methods
        function obj = CenterOfRotation(varargin)
            % Constructor
            obj@mufo.core.UFOTask('correlate-stacks');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setSearchRange(obj, range)
            % Set search range for center detection
            %
            % Input:
            %   range - Maximum pixels to search from image center
            
            if range <= 0
                error('CenterOfRotation:InvalidRange', ...
                    'Search range must be positive');
            end
            
            obj.setParameter('search-range', range);
        end
        
        function setSearchStep(obj, step)
            % Set search step size
            %
            % Input:
            %   step - Step size in pixels (can be < 1 for sub-pixel)
            
            if step <= 0
                error('CenterOfRotation:InvalidStep', ...
                    'Search step must be positive');
            end
            
            obj.setParameter('search-step', step);
        end
        
        function setNumProjections(obj, num)
            % Set number of projections in the dataset
            %
            % Input:
            %   num - Total number of projections
            
            if num < 2 || mod(num, 1) ~= 0
                error('CenterOfRotation:InvalidNumber', ...
                    'Number of projections must be integer >= 2');
            end
            
            obj.setParameter('number', num);
        end
        
        function center = getDetectedCenter(obj)
            % Get the detected center of rotation
            %
            % Output:
            %   center - Detected center position in pixels
            
            center = obj.centerResult;
            if isempty(center)
                warning('CenterOfRotation:NoResult', ...
                    'No center detected yet. Run the task first.');
            end
        end
        
        function chain = createCenterDetectionChain(obj, inputPath, angleStep)
            % Create complete chain for center detection
            %
            % Inputs:
            %   inputPath - Path to projection data
            %   angleStep - Angle step between projections (radians)
            % Output:
            %   chain - UFOChain configured for center detection
            
            chain = mufo.core.UFOChain('Name', 'Center of Rotation Detection');
            
            % Read projections
            readTask = mufo.tasks.Read();
            readTask.setPath(inputPath);
            chain.addTask(readTask, 'Name', 'read-projections');
            
            % Calculate indices for 0° and 180° projections
            numAngles = round(pi / angleStep);  % Number of projections in 180°
            
            % Create two branches for opposite projections
            % Branch 1: Extract projection at 0°
            slice1 = mufo.tasks.Slice();
            slice1.setSlice(0);
            chain.addTask(slice1, 'Input', 'read-projections', ...
                         'Name', 'proj-0deg');
            
            % Branch 2: Extract projection at 180°
            slice2 = mufo.tasks.Slice();
            slice2.setSlice(numAngles);
            chain.addTask(slice2, 'Input', 'read-projections', ...
                         'Name', 'proj-180deg');
            
            % Flip the 180° projection horizontally
            flipTask = mufo.tasks.HorizontalFlip();
            chain.addTask(flipTask, 'Input', 'proj-180deg', ...
                         'Name', 'proj-180deg-flipped');
            
            % Stack projections for correlation
            stackTask = mufo.tasks.Stack();
            stackTask.setNumber(2);
            chain.addTask(stackTask, 'Inputs', {'proj-0deg', 'proj-180deg-flipped'}, ...
                         'Name', 'stacked-projections');
            
            % Add correlation task
            obj.setNumProjections(2);
            chain.addTask(obj, 'Input', 'stacked-projections', ...
                         'Name', 'correlate');
            
            % Write results
            writeTask = mufo.tasks.Write();
            writeTask.setPath('center_detection_result.tif');
            chain.addTask(writeTask, 'Input', 'correlate');
        end
        
        function center = detectFromSinogram(obj, sinogram, angleStep)
            % Detect center from a single sinogram
            %
            % Inputs:
            %   sinogram - 2D sinogram data
            %   angleStep - Angle step between projections (radians)
            % Output:
            %   center - Detected center of rotation
            
            [width, numAngles] = size(sinogram);
            
            % Find indices for opposite projections
            idx1 = 1;
            idx2 = round(numAngles / 2) + 1;  % 180 degrees apart
            
            if idx2 > numAngles
                error('CenterOfRotation:InsufficientData', ...
                    'Sinogram does not cover 180 degrees');
            end
            
            % Extract opposite projections
            proj1 = sinogram(:, idx1);
            proj2 = sinogram(:, idx2);
            
            % Get search parameters
            searchRange = obj.getParameter('search-range');
            searchStep = obj.getParameter('search-step');
            
            if isempty(searchRange)
                searchRange = round(width * 0.1);  % 10% of width
            end
            if isempty(searchStep)
                searchStep = 0.5;
            end
            
            % Image center
            imageCenter = width / 2;
            
            % Search positions
            shifts = -searchRange:searchStep:searchRange;
            correlations = zeros(size(shifts));
            
            % Flip proj2 for correlation
            proj2_flipped = flipud(proj2);
            
            % Test each shift
            for i = 1:length(shifts)
                shift = shifts(i);
                
                % Shift the flipped projection
                if shift > 0
                    proj2_shifted = [zeros(round(shift), 1); ...
                                    proj2_flipped(1:end-round(shift))];
                elseif shift < 0
                    proj2_shifted = [proj2_flipped(-round(shift)+1:end); ...
                                    zeros(-round(shift), 1)];
                else
                    proj2_shifted = proj2_flipped;
                end
                
                % Ensure same size
                if length(proj2_shifted) > length(proj1)
                    proj2_shifted = proj2_shifted(1:length(proj1));
                elseif length(proj2_shifted) < length(proj1)
                    proj2_shifted = [proj2_shifted; zeros(length(proj1)-length(proj2_shifted), 1)];
                end
                
                % Calculate correlation
                correlations(i) = corr(proj1, proj2_shifted);
            end
            
            % Find maximum correlation
            [maxCorr, maxIdx] = max(correlations);
            bestShift = shifts(maxIdx);
            
            % Calculate center
            center = imageCenter + bestShift/2;
            obj.centerResult = center;
            
            % Display results
            fprintf('Center of rotation detection:\n');
            fprintf('  Image width: %d pixels\n', width);
            fprintf('  Search range: ±%.1f pixels\n', searchRange);
            fprintf('  Maximum correlation: %.4f\n', maxCorr);
            fprintf('  Detected center: %.2f pixels\n', center);
            
            % Optional: plot correlation curve
            if nargout == 0
                figure('Name', 'Center of Rotation Detection');
                plot(imageCenter + shifts/2, correlations, 'b-', 'LineWidth', 2);
                hold on;
                plot(center, maxCorr, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
                xlabel('Center Position (pixels)');
                ylabel('Correlation');
                title(sprintf('Center of Rotation: %.2f pixels', center));
                grid on;
                xlim([imageCenter - searchRange, imageCenter + searchRange]);
            end
        end
        
        function info = getCorrelationInfo(obj, projectionDims)
            % Get information about correlation operation
            %
            % Input:
            %   projectionDims - Dimensions of projection data [width, height, angles]
            % Output:
            %   info - Struct with correlation information
            
            info = struct();
            
            % Get parameters
            searchRange = obj.getParameter('search-range');
            searchStep = obj.getParameter('search-step');
            
            if isempty(searchRange)
                searchRange = round(projectionDims(1) * 0.1);
            end
            if isempty(searchStep)
                searchStep = 0.5;
            end
            
            % Calculate correlation info
            info.imageWidth = projectionDims(1);
            info.imageHeight = projectionDims(2);
            info.numProjections = projectionDims(3);
            info.searchRange = searchRange;
            info.searchStep = searchStep;
            info.numTests = length(-searchRange:searchStep:searchRange);
            info.imageCenter = projectionDims(1) / 2;
            info.searchMin = info.imageCenter - searchRange;
            info.searchMax = info.imageCenter + searchRange;
            
            % Memory requirements
            info.memoryPerTest = 2 * projectionDims(1) * projectionDims(2) * 4; % 2 images, float32
            info.totalMemory = info.memoryPerTest * info.numTests;
            
            % Store info
            obj.correlationInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== Center Detection Information ===\n');
                fprintf('Image dimensions: %d x %d\n', info.imageWidth, info.imageHeight);
                fprintf('Number of projections: %d\n', info.numProjections);
                fprintf('Search range: %.1f to %.1f pixels\n', info.searchMin, info.searchMax);
                fprintf('Search step: %.2f pixels\n', info.searchStep);
                fprintf('Number of tests: %d\n', info.numTests);
                fprintf('Memory required: %.2f MB\n', info.totalMemory / 1e6);
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
                
                % Set search parameters if available
                if isfield(reconParams, 'centerSearchRange')
                    obj.setSearchRange(reconParams.centerSearchRange);
                end
                
                if isfield(reconParams, 'centerSearchStep')
                    obj.setSearchStep(reconParams.centerSearchStep);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'searchRange')
                    obj.setSearchRange(config.searchRange);
                end
                
                if isfield(config, 'searchStep')
                    obj.setSearchStep(config.searchStep);
                end
                
                if isfield(config, 'numProjections')
                    obj.setNumProjections(config.numProjections);
                end
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define task parameters
            
            % Number of items in second stream (required for correlate-stacks)
            obj.addParameterDef('number', 'integer', ...
                'description', 'Number of projections to correlate', ...
                'required', true);
            
            % Search range (custom parameter)
            obj.addParameterDef('search-range', 'numeric', ...
                'description', 'Maximum pixels to search from center', ...
                'range', [1, inf], ...
                'default', []);
            
            % Search step (custom parameter)
            obj.addParameterDef('search-step', 'numeric', ...
                'description', 'Step size for search (pixels)', ...
                'range', [0.1, inf], ...
                'default', 0.5);
        end
        
        function valid = validate(obj)
            % Validate parameters
            valid = true;
            
            % Check if number is set
            if isempty(obj.getParameter('number'))
                valid = false;
                warning('CenterOfRotation:NoNumber', ...
                    'Number of projections must be set');
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate center of rotation detection
            
            fprintf('\n=== Center of Rotation Detection Demo ===\n\n');
            
            % Create synthetic data
            angles = linspace(0, pi, 180);
            detector = linspace(-1, 1, 512);
            [D, A] = meshgrid(detector, angles);
            
            % Create phantom with known offset
            trueCenter = 256 + 30;  % 30 pixels off center
            phantom = zeros(512, 180);
            
            for i = 1:180
                theta = angles(i);
                % Circle at (0.15, 0) with radius 0.3
                t = D(i,:) * cos(theta) + 0.15 * sin(theta);
                phantom(:,i) = abs(t) < 0.3;
            end
            
            % Apply center offset
            phantom = circshift(phantom, [30, 0]);
            
            % Add noise
            phantom = phantom + 0.05 * randn(size(phantom));
            
            % Create task and detect center
            task = mufo.tasks.CenterOfRotation();
            task.setSearchRange(50);
            task.setSearchStep(0.25);
            
            % Detect center
            detectedCenter = task.detectFromSinogram(phantom', pi/180);
            
            fprintf('\nTrue center: %.1f pixels\n', trueCenter);
            fprintf('Detected center: %.2f pixels\n', detectedCenter);
            fprintf('Error: %.2f pixels\n', abs(detectedCenter - trueCenter));
            
            % Show sinogram
            figure('Name', 'Test Sinogram');
            imagesc(phantom');
            colormap gray;
            xlabel('Detector Position');
            ylabel('Projection Angle');
            title('Synthetic Sinogram with Off-Center Object');
            hold on;
            plot([trueCenter, trueCenter], [1, 180], 'r--', 'LineWidth', 2);
            plot([detectedCenter, detectedCenter], [1, 180], 'g--', 'LineWidth', 2);
            legend('True Center', 'Detected Center');
            
            fprintf('\nDemo complete.\n');
        end
    end
end