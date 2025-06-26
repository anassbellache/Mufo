classdef FFT < mufo.core.UFOTask
    % FFT - UFO FFT task wrapper
    %
    % Performs Fast Fourier Transform along specified dimensions.
    % Essential component for Fourier-based reconstruction methods.
    %
    % Example:
    %   task = mufo.tasks.FFT();
    %   task.setDimensions(1);        % 1D FFT along first dimension
    %   task.setSize([2048, 2048]);   % Specify FFT size
    
    properties (Access = private)
        fftInfo         % Information about FFT operation
    end
    
    methods
        function obj = FFT(varargin)
            % Constructor
            obj@mufo.core.UFOTask('fft');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setDimensions(obj, dims)
            % Set FFT dimensions (1D, 2D, or 3D)
            %
            % Input:
            %   dims - Number of dimensions (1, 2, or 3)
            
            if ~any(dims == [1, 2, 3])
                error('FFT:InvalidDimensions', ...
                    'Dimensions must be 1, 2, or 3');
            end
            
            obj.setParameter('dimensions', dims);
        end
        
        function setSize(obj, sizes)
            % Set FFT sizes for each dimension
            %
            % Input:
            %   sizes - Array of sizes [size_x] or [size_x, size_y] or [size_x, size_y, size_z]
            
            if length(sizes) >= 1
                obj.setParameter('size-x', sizes(1));
            end
            if length(sizes) >= 2
                obj.setParameter('size-y', sizes(2));
            end
            if length(sizes) >= 3
                obj.setParameter('size-z', sizes(3));
            end
        end
        
        function setPadding(obj, enable)
            % Enable/disable zero padding
            %
            % Input:
            %   enable - true to enable padding, false otherwise
            
            obj.setParameter('padding', enable);
        end
        
        function setInPlace(obj, enable)
            % Enable/disable in-place FFT
            %
            % Input:
            %   enable - true for in-place operation
            
            obj.setParameter('in-place', enable);
        end
        
        function info = getFFTInfo(obj, inputDims)
            % Get information about FFT operation
            %
            % Input:
            %   inputDims - Input data dimensions
            % Output:
            %   info - Structure with FFT information
            
            dims = obj.getParameter('dimensions');
            if isempty(dims)
                dims = 1;  % Default 1D
            end
            
            info = struct();
            info.inputDims = inputDims;
            info.dimensions = dims;
            
            % Get FFT sizes
            info.fftSize = zeros(1, dims);
            if dims >= 1
                sizeX = obj.getParameter('size-x');
                info.fftSize(1) = ifelse(isempty(sizeX), inputDims(1), sizeX);
            end
            if dims >= 2
                sizeY = obj.getParameter('size-y');
                info.fftSize(2) = ifelse(isempty(sizeY), inputDims(2), sizeY);
            end
            if dims >= 3
                sizeZ = obj.getParameter('size-z');
                info.fftSize(3) = ifelse(isempty(sizeZ), inputDims(3), sizeZ);
            end
            
            % Calculate padding
            info.padding = zeros(1, dims);
            for i = 1:dims
                info.padding(i) = info.fftSize(i) - inputDims(i);
            end
            
            % Memory estimate
            info.complexOutput = true;
            info.bytesPerElement = 8;  % Complex single precision
            info.memoryRequired = prod(info.fftSize) * info.bytesPerElement;
            
            obj.fftInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== FFT Operation ===\n');
                fprintf('Dimensions: %dD\n', dims);
                fprintf('Input size: [%s]\n', num2str(inputDims));
                fprintf('FFT size: [%s]\n', num2str(info.fftSize));
                if any(info.padding > 0)
                    fprintf('Padding: [%s]\n', num2str(info.padding));
                end
                fprintf('Memory required: %.2f MB\n', info.memoryRequired / 1e6);
            end
        end
        
        function configureForReconstruction(obj, detectorWidth, numProjections)
            % Configure FFT for typical reconstruction workflow
            %
            % Inputs:
            %   detectorWidth - Width of detector
            %   numProjections - Number of projections
            
            % For FBP, we typically do 1D FFT along detector dimension
            obj.setDimensions(1);
            
            % Set size to next power of 2 for efficiency
            fftSize = 2^nextpow2(detectorWidth);
            obj.setSize(fftSize);
            
            fprintf('Configured FFT for reconstruction:\n');
            fprintf('  Detector width: %d -> FFT size: %d\n', ...
                detectorWidth, fftSize);
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get reconstruction settings
                algorithm = config.get('reconstruction.algorithm', 'fbp');
                
                % Configure based on algorithm
                if strcmpi(algorithm, 'fbp')
                    % FBP uses 1D FFT along detector dimension
                    obj.setDimensions(1);
                end
                
                % Check for FFT-specific settings
                fftConfig = config.get('processing.fft', struct());
                if isfield(fftConfig, 'padding')
                    obj.setPadding(fftConfig.padding);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'dimensions')
                    obj.setDimensions(config.dimensions);
                end
                
                if isfield(config, 'size')
                    obj.setSize(config.size);
                end
                
                if isfield(config, 'padding')
                    obj.setPadding(config.padding);
                end
                
                if isfield(config, 'inPlace')
                    obj.setInPlace(config.inPlace);
                end
            end
        end
        
        function validateForInput(obj, inputDims)
            % Validate FFT configuration for given input
            %
            % Input:
            %   inputDims - Input data dimensions
            
            dims = obj.getParameter('dimensions');
            if isempty(dims)
                dims = 1;
            end
            
            % Check dimension compatibility
            if dims > length(inputDims)
                error('FFT:DimensionMismatch', ...
                    'FFT dimensions (%d) exceed data dimensions (%d)', ...
                    dims, length(inputDims));
            end
            
            % Check sizes
            info = obj.getFFTInfo(inputDims);
            
            % Warn about excessive padding
            for i = 1:dims
                if info.padding(i) > inputDims(i)
                    warning('FFT:ExcessivePadding', ...
                        'Dimension %d: padding (%d) exceeds data size (%d)', ...
                        i, info.padding(i), inputDims(i));
                end
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define FFT parameters
            
            % Number of dimensions
            obj.addParameterDef('dimensions', 'integer', ...
                'description', 'Number of FFT dimensions', ...
                'values', [1, 2, 3], ...
                'default', 1);
            
            % FFT sizes
            obj.addParameterDef('size-x', 'integer', ...
                'description', 'FFT size in X dimension', ...
                'range', [1, inf]);
            
            obj.addParameterDef('size-y', 'integer', ...
                'description', 'FFT size in Y dimension', ...
                'range', [1, inf]);
            
            obj.addParameterDef('size-z', 'integer', ...
                'description', 'FFT size in Z dimension', ...
                'range', [1, inf]);
            
            % Options
            obj.addParameterDef('padding', 'logical', ...
                'description', 'Enable zero padding', ...
                'default', true);
            
            obj.addParameterDef('in-place', 'logical', ...
                'description', 'Perform FFT in-place', ...
                'default', false);
            
            % Normalization
            obj.addParameterDef('normalize', 'logical', ...
                'description', 'Normalize FFT output', ...
                'default', false);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            dims = obj.getParameter('dimensions');
            if ~isempty(dims)
                % Check that appropriate sizes are set
                if dims >= 1 && isempty(obj.getParameter('size-x'))
                    % Size not required - will use input size
                end
                
                % Check for valid size values
                sizes = [obj.getParameter('size-x'), ...
                        obj.getParameter('size-y'), ...
                        obj.getParameter('size-z')];
                
                for i = 1:length(sizes)
                    if ~isempty(sizes(i)) && sizes(i) <= 0
                        valid = false;
                        warning('FFT:InvalidSize', ...
                            'FFT size must be positive');
                        break;
                    end
                end
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
            % Demonstrate FFT operations
            
            fprintf('\n=== FFT Task Demo ===\n\n');
            
            % Create task
            task = mufo.tasks.FFT();
            
            % Example 1: 1D FFT for FBP
            fprintf('1. 1D FFT for FBP reconstruction:\n');
            task.setDimensions(1);
            task.configureForReconstruction(2048, 1800);
            task.getFFTInfo([2048, 1800]);
            
            % Example 2: 2D FFT
            fprintf('\n2. 2D FFT with custom size:\n');
            task.setDimensions(2);
            task.setSize([4096, 4096]);
            task.getFFTInfo([3000, 3000]);
            
            % Example 3: 3D FFT with padding
            fprintf('\n3. 3D FFT with automatic padding:\n');
            task.setDimensions(3);
            task.setSize([]);  % Clear custom sizes
            task.setPadding(true);
            task.getFFTInfo([100, 100, 100]);
            
            fprintf('\nDemo complete.\n');
        end
    end
end