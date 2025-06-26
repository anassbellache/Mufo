classdef RemoveRings < mufo.core.UFOTask
    % RemoveRings - UFO ring removal task wrapper
    %
    % Removes ring artifacts from sinogram data using various methods.
    % Ring artifacts appear as vertical stripes in sinograms.
    %
    % Example:
    %   task = mufo.tasks.RemoveRings();
    %   task.setMethod('wavelet');
    %   task.setThreshold(1.0);
    %   task.setSigma(1.0);
    
    properties (Access = private)
        ringInfo        % Information about ring removal
    end
    
    methods
        function obj = RemoveRings(varargin)
            % Constructor
            obj@mufo.core.UFOTask('remove-rings');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setMethod(obj, method)
            % Set ring removal method
            %
            % Input:
            %   method - 'wavelet' or 'morphological'
            
            validMethods = {'wavelet', 'morphological'};
            if ~any(strcmpi(method, validMethods))
                error('RemoveRings:InvalidMethod', ...
                    'Method must be wavelet or morphological');
            end
            
            obj.setParameter('method', method);
        end
        
        function setThreshold(obj, threshold)
            % Set threshold for ring detection
            %
            % Input:
            %   threshold - Detection threshold (typically 0.5-2.0)
            
            if threshold <= 0
                error('RemoveRings:InvalidThreshold', ...
                    'Threshold must be positive');
            end
            
            obj.setParameter('threshold', threshold);
        end
        
        function setSigma(obj, sigma)
            % Set Gaussian sigma for filtering
            %
            % Input:
            %   sigma - Gaussian sigma (typically 1-5)
            
            if sigma <= 0
                error('RemoveRings:InvalidSigma', ...
                    'Sigma must be positive');
            end
            
            obj.setParameter('sigma', sigma);
        end
        
        function setWaveletLevel(obj, level)
            % Set wavelet decomposition level
            %
            % Input:
            %   level - Decomposition level (1-5)
            
            if level < 1 || level > 5 || mod(level, 1) ~= 0
                error('RemoveRings:InvalidLevel', ...
                    'Level must be integer between 1 and 5');
            end
            
            obj.setParameter('wavelet-level', level);
        end
        
        function setFilterSize(obj, size)
            % Set morphological filter size
            %
            % Input:
            %   size - Filter size (odd integer)
            
            if mod(size, 2) == 0
                error('RemoveRings:InvalidSize', ...
                    'Filter size must be odd');
            end
            
            obj.setParameter('filter-size', size);
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get ring removal settings
                ringConfig = config.get('preprocessing.ringRemoval', struct());
                
                if isfield(ringConfig, 'enabled') && ~ringConfig.enabled
                    warning('RemoveRings:Disabled', ...
                        'Ring removal is disabled in configuration');
                end
                
                if isfield(ringConfig, 'method')
                    obj.setMethod(ringConfig.method);
                end
                
                % Get method-specific parameters
                params = ringConfig.parameters;
                if isstruct(params)
                    if isfield(params, 'threshold')
                        obj.setThreshold(params.threshold);
                    end
                    if isfield(params, 'sigma')
                        obj.setSigma(params.sigma);
                    end
                    if isfield(params, 'waveletLevel')
                        obj.setWaveletLevel(params.waveletLevel);
                    end
                    if isfield(params, 'filterSize')
                        obj.setFilterSize(params.filterSize);
                    end
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'method')
                    obj.setMethod(config.method);
                end
                
                if isfield(config, 'threshold')
                    obj.setThreshold(config.threshold);
                end
                
                if isfield(config, 'sigma')
                    obj.setSigma(config.sigma);
                end
                
                if isfield(config, 'waveletLevel')
                    obj.setWaveletLevel(config.waveletLevel);
                end
                
                if isfield(config, 'filterSize')
                    obj.setFilterSize(config.filterSize);
                end
            end
        end
        
        function analyzeRings(obj, sinogram)
            % Analyze ring artifacts in sinogram
            %
            % Input:
            %   sinogram - 2D sinogram data
            
            [width, numAngles] = size(sinogram);
            
            % Compute mean projection
            meanProj = mean(sinogram, 2);
            
            % Compute standard deviation along angles
            stdProj = std(sinogram, 0, 2);
            
            % Detect potential ring positions
            % Rings appear as consistent features across all angles
            threshold = obj.getParameter('threshold');
            if isempty(threshold)
                threshold = 1.0;
            end
            
            % Normalize standard deviation
            stdNorm = stdProj / mean(stdProj);
            ringPositions = find(stdNorm < (1/threshold));
            
            % Create figure
            figure('Name', 'Ring Artifact Analysis');
            
            % Show sinogram
            subplot(2,2,1);
            imagesc(sinogram);
            colormap gray;
            xlabel('Projection Angle');
            ylabel('Detector Position');
            title('Original Sinogram');
            colorbar;
            
            % Show mean projection
            subplot(2,2,2);
            plot(meanProj, 1:width, 'b-');
            hold on;
            plot(meanProj(ringPositions), ringPositions, 'ro');
            xlabel('Mean Intensity');
            ylabel('Detector Position');
            title('Mean Projection');
            grid on;
            
            % Show standard deviation
            subplot(2,2,3);
            plot(stdProj, 1:width, 'b-');
            hold on;
            plot(stdProj(ringPositions), ringPositions, 'ro');
            xlabel('Standard Deviation');
            ylabel('Detector Position');
            title('Intensity Variation');
            grid on;
            
            % Show normalized std
            subplot(2,2,4);
            plot(stdNorm, 1:width, 'b-');
            hold on;
            plot([0, width], [1/threshold, 1/threshold], 'r--');
            xlabel('Normalized Std Dev');
            ylabel('Detector Position');
            title('Ring Detection');
            legend('Std Dev', sprintf('Threshold (1/%.1f)', threshold));
            grid on;
            
            % Store ring info
            obj.ringInfo = struct(...
                'numRings', length(ringPositions), ...
                'positions', ringPositions, ...
                'severity', 1 - stdNorm(ringPositions));
            
            fprintf('\nRing Analysis Results:\n');
            fprintf('  Detected rings: %d\n', obj.ringInfo.numRings);
            if obj.ringInfo.numRings > 0
                fprintf('  Ring positions: %s\n', ...
                    num2str(ringPositions(1:min(10, end))'));
                if obj.ringInfo.numRings > 10
                    fprintf('  ... and %d more\n', obj.ringInfo.numRings - 10);
                end
            end
        end
        
        function previewRemoval(obj, sinogram)
            % Preview ring removal on sample sinogram
            %
            % Input:
            %   sinogram - 2D sinogram data
            
            % This is a simplified preview - actual UFO implementation
            % will be more sophisticated
            
            method = obj.getParameter('method');
            if isempty(method)
                method = 'wavelet';
            end
            
            fprintf('Simulating %s ring removal...\n', method);
            
            % Simple simulation based on method
            switch lower(method)
                case 'wavelet'
                    % Simple high-pass filtering simulation
                    sigma = obj.getParameter('sigma');
                    if isempty(sigma), sigma = 1.0; end
                    
                    % Apply Gaussian filter along angle dimension
                    filtered = sinogram;
                    for i = 1:size(sinogram, 1)
                        filtered(i, :) = sinogram(i, :) - ...
                            imgaussfilt(sinogram(i, :), sigma);
                    end
                    corrected = sinogram - filtered;
                    
                case 'morphological'
                    % Simple morphological filtering simulation
                    filterSize = obj.getParameter('filter-size');
                    if isempty(filterSize), filterSize = 5; end
                    
                    % Apply opening along detector dimension
                    se = strel('line', filterSize, 90);
                    corrected = imopen(sinogram, se);
                    
                otherwise
                    corrected = sinogram;
            end
            
            % Display results
            figure('Name', 'Ring Removal Preview');
            
            subplot(1,3,1);
            imagesc(sinogram);
            colormap gray;
            title('Original');
            xlabel('Angle');
            ylabel('Detector');
            colorbar;
            
            subplot(1,3,2);
            imagesc(corrected);
            colormap gray;
            title(sprintf('After %s', method));
            xlabel('Angle');
            ylabel('Detector');
            colorbar;
            
            subplot(1,3,3);
            imagesc(sinogram - corrected);
            colormap gray;
            title('Removed Features');
            xlabel('Angle');
            ylabel('Detector');
            colorbar;
            
            % Show profiles
            figure('Name', 'Ring Removal Profiles');
            midAngle = round(size(sinogram, 2) / 2);
            
            plot(sinogram(:, midAngle), 'b-', 'LineWidth', 1.5);
            hold on;
            plot(corrected(:, midAngle), 'r-', 'LineWidth', 1.5);
            xlabel('Detector Position');
            ylabel('Intensity');
            title(sprintf('Profile at Angle %d', midAngle));
            legend('Original', 'Corrected');
            grid on;
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define ring removal parameters
            
            % Method
            obj.addParameterDef('method', 'string', ...
                'description', 'Ring removal method', ...
                'values', {'wavelet', 'morphological'}, ...
                'default', 'wavelet');
            
            % Detection threshold
            obj.addParameterDef('threshold', 'numeric', ...
                'description', 'Detection threshold', ...
                'range', [0, inf], ...
                'default', 1.0);
            
            % Gaussian sigma
            obj.addParameterDef('sigma', 'numeric', ...
                'description', 'Gaussian sigma for filtering', ...
                'range', [0, inf], ...
                'default', 1.0);
            
            % Wavelet parameters
            obj.addParameterDef('wavelet-level', 'integer', ...
                'description', 'Wavelet decomposition level', ...
                'range', [1, 5], ...
                'default', 3);
            
            obj.addParameterDef('wavelet-type', 'string', ...
                'description', 'Wavelet type', ...
                'values', {'db1', 'db2', 'db3', 'db4', 'sym4', 'sym8'}, ...
                'default', 'db3');
            
            % Morphological parameters
            obj.addParameterDef('filter-size', 'integer', ...
                'description', 'Morphological filter size', ...
                'range', [3, 21], ...
                'default', 5);
            
            % Processing options
            obj.addParameterDef('preserve-edges', 'logical', ...
                'description', 'Preserve edge information', ...
                'default', true);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            method = obj.getParameter('method');
            
            if strcmpi(method, 'morphological')
                filterSize = obj.getParameter('filter-size');
                if ~isempty(filterSize) && mod(filterSize, 2) == 0
                    valid = false;
                    warning('RemoveRings:EvenFilterSize', ...
                        'Morphological filter size must be odd');
                end
            end
            
            % Check for reasonable threshold
            threshold = obj.getParameter('threshold');
            if ~isempty(threshold)
                if threshold < 0.1
                    warning('RemoveRings:LowThreshold', ...
                        'Very low threshold may remove valid features');
                elseif threshold > 5
                    warning('RemoveRings:HighThreshold', ...
                        'Very high threshold may miss ring artifacts');
                end
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate ring removal operations
            
            fprintf('\n=== Ring Removal Task Demo ===\n\n');
            
            % Create synthetic sinogram with rings
            angles = linspace(0, 2*pi, 360);
            detectorSize = 512;
            
            % Create base sinogram
            sinogram = zeros(detectorSize, length(angles));
            for i = 1:length(angles)
                % Simple phantom projection
                t = linspace(-1, 1, detectorSize);
                sinogram(:, i) = exp(-t.^2 / 0.1);
            end
            
            % Add ring artifacts
            ringPositions = [100, 200, 300, 400];
            ringIntensities = [0.1, -0.15, 0.2, -0.1];
            
            for i = 1:length(ringPositions)
                sinogram(ringPositions(i), :) = ...
                    sinogram(ringPositions(i), :) + ringIntensities(i);
            end
            
            % Add some noise
            sinogram = sinogram + 0.02 * randn(size(sinogram));
            
            % Create task
            task = mufo.tasks.RemoveRings();
            
            % Example 1: Analyze rings
            fprintf('1. Analyzing ring artifacts:\n');
            task.setThreshold(1.5);
            task.analyzeRings(sinogram);
            
            % Example 2: Wavelet method
            fprintf('\n2. Wavelet ring removal:\n');
            task.setMethod('wavelet');
            task.setSigma(2.0);
            task.previewRemoval(sinogram);
            
            % Example 3: Morphological method
            fprintf('\n3. Morphological ring removal:\n');
            task.setMethod('morphological');
            task.setFilterSize(7);
            task.previewRemoval(sinogram);
            
            fprintf('\nDemo complete.\n');
        end
    end
end