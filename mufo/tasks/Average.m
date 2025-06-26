classdef Average < mufo.core.UFOTask
    % Average - UFO averaging task wrapper
    %
    % Averages multiple frames together. Commonly used for creating
    % averaged dark and flat field references from multiple acquisitions.
    %
    % Example:
    %   task = mufo.tasks.Average();
    %   task.setMethod('mean');  % or 'median'
    
    properties (Access = private)
        averageInfo     % Information about averaging operation
    end
    
    methods
        function obj = Average(varargin)
            % Constructor
            obj@mufo.core.UFOTask('average');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setMethod(obj, method)
            % Set averaging method
            %
            % Input:
            %   method - 'mean' or 'median'
            
            validMethods = {'mean', 'median'};
            if ~any(strcmpi(method, validMethods))
                error('Average:InvalidMethod', ...
                    'Method must be mean or median');
            end
            
            obj.setParameter('method', method);
        end
        
        function setNumInputs(obj, num)
            % Set expected number of inputs
            %
            % Input:
            %   num - Number of frames to average
            
            if num < 1 || mod(num, 1) ~= 0
                error('Average:InvalidNumber', ...
                    'Number of inputs must be positive integer');
            end
            
            obj.setParameter('num-inputs', num);
        end
        
        function setWeights(obj, weights)
            % Set weights for weighted averaging
            %
            % Input:
            %   weights - Array of weights (must sum to 1)
            
            if any(weights < 0)
                error('Average:InvalidWeights', ...
                    'Weights must be non-negative');
            end
            
            % Normalize weights
            weights = weights / sum(weights);
            
            % Convert to string for UFO
            weightStr = sprintf('%.6f,', weights);
            weightStr(end) = [];  % Remove trailing comma
            
            obj.setParameter('weights', weightStr);
        end
        
        function setRobust(obj, enable)
            % Enable robust averaging (outlier rejection)
            %
            % Input:
            %   enable - true for robust averaging
            
            obj.setParameter('robust', enable);
        end
        
        function info = getAverageInfo(obj, inputDims, numFrames)
            % Get information about averaging operation
            %
            % Inputs:
            %   inputDims - Dimensions of single frame
            %   numFrames - Number of frames to average
            % Output:
            %   info - Structure with averaging information
            
            info = struct();
            info.inputDims = inputDims;
            info.numFrames = numFrames;
            info.outputDims = inputDims;  % Same as input
            
            % Get method
            method = obj.getParameter('method');
            if isempty(method)
                method = 'mean';
            end
            info.method = method;
            
            % Get weights if specified
            weightsStr = obj.getParameter('weights');
            if ~isempty(weightsStr)
                info.weights = sscanf(weightsStr, '%f,');
                info.weighted = true;
            else
                info.weights = ones(numFrames, 1) / numFrames;
                info.weighted = false;
            end
            
            % Memory estimate
            if strcmpi(method, 'median')
                % Median requires all frames in memory
                info.memoryRequired = prod(inputDims) * numFrames * 4;
            else
                % Mean can be computed incrementally
                info.memoryRequired = prod(inputDims) * 4 * 2;  % Input + accumulator
            end
            
            obj.averageInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== Averaging Operation ===\n');
                fprintf('Method: %s\n', method);
                fprintf('Input frames: %d × [%s]\n', numFrames, num2str(inputDims));
                fprintf('Output: [%s]\n', num2str(info.outputDims));
                if info.weighted
                    fprintf('Weighted averaging enabled\n');
                end
                fprintf('Memory required: %.2f MB\n', info.memoryRequired / 1e6);
            end
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Check for averaging settings in preprocessing
                avgConfig = config.get('preprocessing.averaging', struct());
                
                if isfield(avgConfig, 'method')
                    obj.setMethod(avgConfig.method);
                end
                
                if isfield(avgConfig, 'robust')
                    obj.setRobust(avgConfig.robust);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'method')
                    obj.setMethod(config.method);
                end
                
                if isfield(config, 'numInputs')
                    obj.setNumInputs(config.numInputs);
                end
                
                if isfield(config, 'weights')
                    obj.setWeights(config.weights);
                end
                
                if isfield(config, 'robust')
                    obj.setRobust(config.robust);
                end
            end
        end
        
        function chain = createDarkFlatAveraging(obj, darkPath, flatPath)
            % Create chain for dark/flat field averaging
            %
            % Inputs:
            %   darkPath - Path to dark fields
            %   flatPath - Path to flat fields
            % Output:
            %   chain - UFOChain with parallel averaging
            
            chain = mufo.core.UFOChain('Name', 'Dark/Flat Averaging');
            
            % Read tasks
            darkRead = mufo.tasks.Read();
            darkRead.setPath(darkPath);
            
            flatRead = mufo.tasks.Read();
            flatRead.setPath(flatPath);
            
            % Average tasks
            darkAvg = mufo.tasks.Average();
            flatAvg = mufo.tasks.Average();
            
            % Build parallel chains
            darkChain = {darkRead, darkAvg};
            flatChain = {flatRead, flatAvg};
            
            chain.addParallel({darkChain, flatChain});
        end
        
        function demonstrateAveraging(obj, data)
            % Demonstrate different averaging methods
            %
            % Input:
            %   data - 3D data array (x, y, frames)
            
            [height, width, numFrames] = size(data);
            
            % Compute different averages
            meanResult = mean(data, 3);
            medianResult = median(data, 3);
            
            % Robust mean (trim 10% outliers)
            sortedData = sort(data, 3);
            trimFrames = round(0.1 * numFrames);
            robustData = sortedData(:, :, trimFrames+1:end-trimFrames);
            robustMean = mean(robustData, 3);
            
            % Display results
            figure('Name', 'Averaging Methods Comparison');
            
            subplot(2,3,1);
            imagesc(data(:,:,1));
            colormap gray;
            title('Single Frame');
            colorbar;
            
            subplot(2,3,2);
            imagesc(meanResult);
            colormap gray;
            title(sprintf('Mean (%d frames)', numFrames));
            colorbar;
            
            subplot(2,3,3);
            imagesc(medianResult);
            colormap gray;
            title('Median');
            colorbar;
            
            subplot(2,3,4);
            imagesc(robustMean);
            colormap gray;
            title('Robust Mean (10% trim)');
            colorbar;
            
            % Show noise reduction
            subplot(2,3,5);
            noiseOriginal = std(data, 0, 3);
            imagesc(noiseOriginal);
            colormap hot;
            title('Original Noise (Std Dev)');
            colorbar;
            
            subplot(2,3,6);
            % Estimate noise reduction
            centerROI = data(height/2-50:height/2+50, width/2-50:width/2+50, :);
            originalNoise = std(centerROI(:));
            meanNoise = std(meanResult(height/2-50:height/2+50, width/2-50:width/2+50));
            
            text(0.1, 0.8, sprintf('Original noise: %.3f', originalNoise));
            text(0.1, 0.6, sprintf('After averaging: %.3f', meanNoise));
            text(0.1, 0.4, sprintf('Noise reduction: %.1f×', ...
                originalNoise / meanNoise));
            text(0.1, 0.2, sprintf('Expected: %.1f×', sqrt(numFrames)));
            axis off;
            title('Noise Statistics');
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define averaging parameters
            
            % Method
            obj.addParameterDef('method', 'string', ...
                'description', 'Averaging method', ...
                'values', {'mean', 'median'}, ...
                'default', 'mean');
            
            % Number of inputs
            obj.addParameterDef('num-inputs', 'integer', ...
                'description', 'Number of inputs to average', ...
                'range', [1, inf]);
            
            % Weights for weighted averaging
            obj.addParameterDef('weights', 'string', ...
                'description', 'Comma-separated weights');
            
            % Robust averaging
            obj.addParameterDef('robust', 'logical', ...
                'description', 'Enable robust averaging', ...
                'default', false);
            
            % Outlier threshold (for robust averaging)
            obj.addParameterDef('outlier-threshold', 'numeric', ...
                'description', 'Outlier rejection threshold', ...
                'range', [0, inf], ...
                'default', 3.0);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            % Check weights if specified
            weightsStr = obj.getParameter('weights');
            if ~isempty(weightsStr)
                try
                    weights = sscanf(weightsStr, '%f,');
                    if any(weights < 0)
                        valid = false;
                        warning('Average:NegativeWeights', ...
                            'Weights must be non-negative');
                    end
                    
                    % Check if weights sum to approximately 1
                    sumWeights = sum(weights);
                    if abs(sumWeights - 1) > 1e-6
                        warning('Average:WeightSum', ...
                            'Weights sum to %.6f (should be 1)', sumWeights);
                    end
                catch
                    valid = false;
                    warning('Average:InvalidWeightFormat', ...
                        'Cannot parse weights string');
                end
            end
            
            % Check if number of inputs matches weights
            numInputs = obj.getParameter('num-inputs');
            if ~isempty(numInputs) && ~isempty(weightsStr)
                weights = sscanf(weightsStr, '%f,');
                if length(weights) ~= numInputs
                    valid = false;
                    warning('Average:WeightMismatch', ...
                        'Number of weights (%d) does not match inputs (%d)', ...
                        length(weights), numInputs);
                end
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate averaging operations
            
            fprintf('\n=== Average Task Demo ===\n\n');
            
            % Create synthetic data with noise
            baseImage = phantom(256);
            noise_level = 0.1;
            numFrames = 10;
            
            data = zeros(256, 256, numFrames);
            for i = 1:numFrames
                data(:,:,i) = baseImage + noise_level * randn(256, 256);
            end
            
            % Create task
            task = mufo.tasks.Average();
            
            % Example 1: Basic averaging info
            fprintf('1. Basic averaging operation:\n');
            task.setMethod('mean');
            task.getAverageInfo([256, 256], numFrames);
            
            % Example 2: Weighted averaging
            fprintf('\n2. Weighted averaging:\n');
            weights = exp(-0.5 * (0:numFrames-1));  % Exponential weights
            weights = weights / sum(weights);
            task.setWeights(weights);
            task.getAverageInfo([256, 256], numFrames);
            
            % Example 3: Compare methods
            fprintf('\n3. Comparing averaging methods:\n');
            task.demonstrateAveraging(data);
            
            % Example 4: Create dark/flat chain
            fprintf('\n4. Creating dark/flat averaging chain:\n');
            chain = task.createDarkFlatAveraging('darks/*.tif', 'flats/*.tif');
            fprintf('Chain created with parallel dark/flat averaging\n');
            
            fprintf('\nDemo complete.\n');
        end
    end
end