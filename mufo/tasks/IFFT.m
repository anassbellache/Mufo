classdef IFFT < mufo.core.UFOTask
    % IFFT - UFO inverse FFT task wrapper
    %
    % Performs Inverse Fast Fourier Transform along specified dimensions.
    % Typically used after filtering in Fourier space for FBP reconstruction.
    %
    % Example:
    %   task = mufo.tasks.IFFT();
    %   task.setDimensions(1);      % 1D IFFT
    %   task.setCrop(2048);         % Crop to original size
    
    properties (Access = private)
        ifftInfo        % Information about IFFT operation
    end
    
    methods
        function obj = IFFT(varargin)
            % Constructor
            obj@mufo.core.UFOTask('ifft');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setDimensions(obj, dims)
            % Set IFFT dimensions (1D, 2D, or 3D)
            %
            % Input:
            %   dims - Number of dimensions (1, 2, or 3)
            
            if ~any(dims == [1, 2, 3])
                error('IFFT:InvalidDimensions', ...
                    'Dimensions must be 1, 2, or 3');
            end
            
            obj.setParameter('dimensions', dims);
        end
        
        function setCrop(obj, sizes)
            % Set output size (crop after IFFT)
            %
            % Input:
            %   sizes - Scalar or array of sizes [size_x] or [size_x, size_y] or [size_x, size_y, size_z]
            
            if isscalar(sizes)
                obj.setParameter('crop-width', sizes);
            else
                if length(sizes) >= 1
                    obj.setParameter('crop-width', sizes(1));
                end
                if length(sizes) >= 2
                    obj.setParameter('crop-height', sizes(2));
                end
                if length(sizes) >= 3
                    obj.setParameter('crop-depth', sizes(3));
                end
            end
        end
        
        function setInPlace(obj, enable)
            % Enable/disable in-place IFFT
            %
            % Input:
            %   enable - true for in-place operation
            
            obj.setParameter('in-place', enable);
        end
        
        function setNormalize(obj, enable)
            % Enable/disable normalization
            %
            % Input:
            %   enable - true to normalize output
            
            obj.setParameter('normalize', enable);
        end
        
        function info = getIFFTInfo(obj, inputDims, fftDims)
            % Get information about IFFT operation
            %
            % Inputs:
            %   inputDims - Input (FFT) data dimensions
            %   fftDims - Original data dimensions before FFT (optional)
            % Output:
            %   info - Structure with IFFT information
            
            dims = obj.getParameter('dimensions');
            if isempty(dims)
                dims = 1;  % Default 1D
            end
            
            info = struct();
            info.inputDims = inputDims;
            info.dimensions = dims;
            
            % Get crop sizes
            cropWidth = obj.getParameter('crop-width');
            cropHeight = obj.getParameter('crop-height');
            cropDepth = obj.getParameter('crop-depth');
            
            % Determine output dimensions
            info.outputDims = inputDims;
            if ~isempty(cropWidth) && dims >= 1
                info.outputDims(1) = cropWidth;
            end
            if ~isempty(cropHeight) && dims >= 2
                info.outputDims(2) = cropHeight;
            end
            if ~isempty(cropDepth) && dims >= 3
                info.outputDims(3) = cropDepth;
            end
            
            % Calculate cropping info
            info.cropRequired = ~isequal(info.inputDims, info.outputDims);
            if info.cropRequired
                info.cropAmount = info.inputDims - info.outputDims;
            else
                info.cropAmount = zeros(size(info.inputDims));
            end
            
            % Memory estimate
            info.realOutput = true;
            info.bytesPerElement = 4;  % Real single precision
            info.memoryRequired = prod(info.outputDims) * info.bytesPerElement;
            
            obj.ifftInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== IFFT Operation ===\n');
                fprintf('Dimensions: %dD\n', dims);
                fprintf('Input size: [%s]\n', num2str(inputDims));
                fprintf('Output size: [%s]\n', num2str(info.outputDims));
                if info.cropRequired
                    fprintf('Cropping: [%s]\n', num2str(info.cropAmount));
                end
                fprintf('Memory required: %.2f MB\n', info.memoryRequired / 1e6);
            end
        end
        
        function configureForReconstruction(obj, fftSize, originalSize)
            % Configure IFFT for typical reconstruction workflow
            %
            % Inputs:
            %   fftSize - Size of FFT output
            %   originalSize - Original data size to crop back to
            
            % For FBP, we typically do 1D IFFT
            obj.setDimensions(1);
            
            % Set cropping to original size
            if fftSize > originalSize
                obj.setCrop(originalSize);
                fprintf('Configured IFFT for reconstruction:\n');
                fprintf('  FFT size: %d -> Original size: %d\n', ...
                    fftSize, originalSize);
            end
            
            % Enable normalization for proper scaling
            obj.setNormalize(true);
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
                    % FBP uses 1D IFFT
                    obj.setDimensions(1);
                    obj.setNormalize(true);
                end
                
                % Check for IFFT-specific settings
                ifftConfig = config.get('processing.ifft', struct());
                if isfield(ifftConfig, 'normalize')
                    obj.setNormalize(ifftConfig.normalize);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'dimensions')
                    obj.setDimensions(config.dimensions);
                end
                
                if isfield(config, 'crop')
                    obj.setCrop(config.crop);
                end
                
                if isfield(config, 'inPlace')
                    obj.setInPlace(config.inPlace);
                end
                
                if isfield(config, 'normalize')
                    obj.setNormalize(config.normalize);
                end
            end
        end
        
        function chain = createFFTIFFTChain(obj, fftTask)
            % Create a complete FFT-IFFT chain for testing
            %
            % Input:
            %   fftTask - Configured FFT task
            % Output:
            %   chain - UFOChain with FFT and IFFT
            
            chain = mufo.core.UFOChain('Name', 'FFT-IFFT Test');
            
            % Add FFT
            chain.addTask(fftTask);
            
            % Configure IFFT to match
            dims = fftTask.getParameter('dimensions');
            obj.setDimensions(dims);
            
            % Add IFFT
            chain.addTask(obj);
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define IFFT parameters
            
            % Number of dimensions
            obj.addParameterDef('dimensions', 'integer', ...
                'description', 'Number of IFFT dimensions', ...
                'values', [1, 2, 3], ...
                'default', 1);
            
            % Crop sizes (to remove padding)
            obj.addParameterDef('crop-width', 'integer', ...
                'description', 'Width after IFFT', ...
                'range', [1, inf]);
            
            obj.addParameterDef('crop-height', 'integer', ...
                'description', 'Height after IFFT', ...
                'range', [1, inf]);
            
            obj.addParameterDef('crop-depth', 'integer', ...
                'description', 'Depth after IFFT', ...
                'range', [1, inf]);
            
            % Options
            obj.addParameterDef('in-place', 'logical', ...
                'description', 'Perform IFFT in-place', ...
                'default', false);
            
            % Normalization
            obj.addParameterDef('normalize', 'logical', ...
                'description', 'Normalize IFFT output', ...
                'default', true);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            dims = obj.getParameter('dimensions');
            
            % Check crop sizes
            cropSizes = [obj.getParameter('crop-width'), ...
                        obj.getParameter('crop-height'), ...
                        obj.getParameter('crop-depth')];
            
            for i = 1:length(cropSizes)
                if ~isempty(cropSizes(i)) && cropSizes(i) <= 0
                    valid = false;
                    warning('IFFT:InvalidCropSize', ...
                        'Crop size must be positive');
                    break;
                end
            end
            
            % Warn if crop dimensions don't match IFFT dimensions
            if ~isempty(dims)
                numCropDims = sum(~cellfun(@isempty, ...
                    {obj.getParameter('crop-width'), ...
                     obj.getParameter('crop-height'), ...
                     obj.getParameter('crop-depth')}));
                
                if numCropDims > 0 && numCropDims < dims
                    warning('IFFT:IncompleteCrop', ...
                        'Crop dimensions (%d) less than IFFT dimensions (%d)', ...
                        numCropDims, dims);
                end
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate IFFT operations
            
            fprintf('\n=== IFFT Task Demo ===\n\n');
            
            % Create tasks
            fftTask = mufo.tasks.FFT();
            ifftTask = mufo.tasks.IFFT();
            
            % Example 1: 1D IFFT for FBP
            fprintf('1. 1D IFFT for FBP reconstruction:\n');
            
            % Configure FFT
            fftTask.setDimensions(1);
            fftTask.setSize(4096);  % Padded size
            
            % Configure matching IFFT
            ifftTask.configureForReconstruction(4096, 2048);
            ifftTask.getIFFTInfo([4096, 1800]);
            
            % Example 2: 2D IFFT with cropping
            fprintf('\n2. 2D IFFT with cropping:\n');
            ifftTask.setDimensions(2);
            ifftTask.setCrop([512, 512]);
            ifftTask.getIFFTInfo([1024, 1024]);
            
            % Example 3: Complete FFT-IFFT chain
            fprintf('\n3. Complete FFT-IFFT chain:\n');
            chain = ifftTask.createFFTIFFTChain(fftTask);
            fprintf('Chain created with %d tasks\n', chain.length);
            
            fprintf('\nDemo complete.\n');
        end
    end
end