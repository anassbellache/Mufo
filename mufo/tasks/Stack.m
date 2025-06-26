classdef Stack < mufo.core.UFOTask
    % Stack - UFO stack task wrapper
    %
    % Stacks multiple 2D inputs into a 3D volume. Used to combine
    % individual slices or projections into a single dataset.
    %
    % Example:
    %   task = mufo.tasks.Stack();
    %   task.setNumber(100);        % Stack 100 slices
    %   task.setDimension(3);       % Stack along Z dimension
    
    properties (Access = private)
        stackInfo       % Information about stacking operation
    end
    
    methods
        function obj = Stack(varargin)
            % Constructor
            obj@mufo.core.UFOTask('stack');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setNumber(obj, num)
            % Set number of items to stack
            %
            % Input:
            %   num - Number of 2D slices to stack
            
            if num < 1 || mod(num, 1) ~= 0
                error('Stack:InvalidNumber', ...
                    'Number must be positive integer');
            end
            
            obj.setParameter('number', num);
        end
        
        function setDimension(obj, dim)
            % Set stacking dimension
            %
            % Input:
            %   dim - Dimension to stack along (1, 2, or 3)
            
            if ~any(dim == [1, 2, 3])
                error('Stack:InvalidDimension', ...
                    'Dimension must be 1, 2, or 3');
            end
            
            obj.setParameter('dimension', dim);
        end
        
        function setBlockSize(obj, size)
            % Set block size for memory-efficient stacking
            %
            % Input:
            %   size - Number of slices per block
            
            if size < 1 || mod(size, 1) ~= 0
                error('Stack:InvalidBlockSize', ...
                    'Block size must be positive integer');
            end
            
            obj.setParameter('block-size', size);
        end
        
        function info = getStackInfo(obj, sliceDims, numSlices)
            % Get information about stacking operation
            %
            % Inputs:
            %   sliceDims - Dimensions of single slice [x, y]
            %   numSlices - Number of slices to stack
            % Output:
            %   info - Structure with stacking information
            
            if nargin < 3
                numSlices = obj.getParameter('number');
                if isempty(numSlices)
                    numSlices = 1;
                end
            end
            
            info = struct();
            info.sliceDims = sliceDims;
            info.numSlices = numSlices;
            
            % Get stacking dimension
            dim = obj.getParameter('dimension');
            if isempty(dim)
                dim = 3;  % Default to Z
            end
            info.dimension = dim;
            
            % Calculate output dimensions
            switch dim
                case 1
                    info.outputDims = [numSlices, sliceDims(1), sliceDims(2)];
                case 2
                    info.outputDims = [sliceDims(1), numSlices, sliceDims(2)];
                case 3
                    info.outputDims = [sliceDims(1), sliceDims(2), numSlices];
            end
            
            % Memory calculation
            info.bytesPerElement = 4;  % Single precision
            info.sliceMemory = prod(sliceDims) * info.bytesPerElement;
            info.totalMemory = info.sliceMemory * numSlices;
            
            % Block processing info
            blockSize = obj.getParameter('block-size');
            if isempty(blockSize)
                blockSize = numSlices;  % Process all at once
            end
            info.blockSize = blockSize;
            info.numBlocks = ceil(numSlices / blockSize);
            info.blockMemory = info.sliceMemory * blockSize;
            
            obj.stackInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== Stack Operation ===\n');
                fprintf('Input: %d slices of [%s]\n', numSlices, num2str(sliceDims));
                fprintf('Output: [%s]\n', num2str(info.outputDims));
                fprintf('Stack dimension: %d\n', dim);
                fprintf('Total memory: %.2f MB\n', info.totalMemory / 1e6);
                if info.numBlocks > 1
                    fprintf('Block processing: %d blocks of %d slices\n', ...
                        info.numBlocks, blockSize);
                    fprintf('Memory per block: %.2f MB\n', info.blockMemory / 1e6);
                end
            end
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get processing configuration
                chunkSize = config.get('processing.chunkSize', []);
                if ~isempty(chunkSize)
                    obj.setBlockSize(chunkSize);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'number')
                    obj.setNumber(config.number);
                end
                
                if isfield(config, 'dimension')
                    obj.setDimension(config.dimension);
                end
                
                if isfield(config, 'blockSize')
                    obj.setBlockSize(config.blockSize);
                end
            end
        end
        
        function chain = createStackingChain(obj, inputPattern, outputFile)
            % Create complete stacking chain
            %
            % Inputs:
            %   inputPattern - File pattern for input slices
            %   outputFile - Output filename
            % Output:
            %   chain - UFOChain for stacking
            
            chain = mufo.core.UFOChain('Name', 'Stacking Operation');
            
            % Read task
            readTask = mufo.tasks.Read();
            readTask.setPath(inputPattern);
            
            % Add to chain
            chain.addTask(readTask);
            chain.addTask(obj);
            
            % Write task
            writeTask = mufo.tasks.Write();
            writeTask.setFilename(outputFile);
            
            % Configure for 3D output
            if endsWith(outputFile, '.h5') || endsWith(outputFile, '.hdf5')
                writeTask.setHDF5Options('/volume/data');
            end
            
            chain.addTask(writeTask);
        end
        
        function demonstrateStacking(obj, data)
            % Demonstrate stacking operations
            %
            % Input:
            %   data - 3D array to demonstrate unstacking/restacking
            
            [height, width, depth] = size(data);
            
            figure('Name', 'Stacking Demonstration');
            
            % Show original volume
            subplot(2,3,1);
            slice = data(:, :, round(depth/2));
            imagesc(slice);
            colormap gray;
            title(sprintf('Original\nMiddle slice'));
            xlabel(sprintf('Volume: %dx%dx%d', height, width, depth));
            
            % Simulate different stacking dimensions
            % Note: This is visualization only - actual UFO operation
            % handles the memory-efficient stacking
            
            % Stack along dimension 1
            subplot(2,3,2);
            reshaped1 = permute(data, [3, 1, 2]);
            imagesc(squeeze(reshaped1(:, :, round(width/2))));
            title('Stack dim 1');
            xlabel(sprintf('Shape: %dx%dx%d', size(reshaped1)));
            
            % Stack along dimension 2
            subplot(2,3,3);
            reshaped2 = permute(data, [1, 3, 2]);
            imagesc(squeeze(reshaped2(:, :, round(width/2))));
            title('Stack dim 2');
            xlabel(sprintf('Shape: %dx%dx%d', size(reshaped2)));
            
            % Show memory usage
            subplot(2,3,4);
            memUsage = [
                height * width * 4;           % Single slice
                height * width * 10 * 4;      % 10 slices
                height * width * depth * 4;   % Full volume
            ] / 1e6;  % Convert to MB
            
            bar(memUsage);
            set(gca, 'XTickLabel', {'1 slice', '10 slices', 'Full volume'});
            ylabel('Memory (MB)');
            title('Memory Requirements');
            grid on;
            
            % Show block processing
            subplot(2,3,5);
            blockSizes = [10, 50, 100, depth];
            numBlocks = ceil(depth ./ blockSizes);
            
            plot(blockSizes, numBlocks, 'bo-', 'LineWidth', 2);
            xlabel('Block Size');
            ylabel('Number of Blocks');
            title('Block Processing');
            grid on;
            
            % Info text
            subplot(2,3,6);
            axis off;
            text(0.1, 0.8, 'Stacking Benefits:');
            text(0.1, 0.6, '• Combine individual files');
            text(0.1, 0.5, '• Memory-efficient processing');
            text(0.1, 0.4, '• Flexible dimension ordering');
            text(0.1, 0.2, 'Block processing allows');
            text(0.1, 0.1, 'handling large datasets');
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define stack parameters
            
            % Number of items
            obj.addParameterDef('number', 'integer', ...
                'description', 'Number of items to stack', ...
                'range', [1, inf], ...
                'default', 1);
            
            % Stacking dimension
            obj.addParameterDef('dimension', 'integer', ...
                'description', 'Dimension to stack along', ...
                'values', [1, 2, 3], ...
                'default', 3);
            
            % Block size for memory efficiency
            obj.addParameterDef('block-size', 'integer', ...
                'description', 'Number of items per processing block', ...
                'range', [1, inf]);
            
            % Options
            obj.addParameterDef('check-dimensions', 'logical', ...
                'description', 'Verify all slices have same dimensions', ...
                'default', true);
            
            obj.addParameterDef('fill-value', 'numeric', ...
                'description', 'Value for missing slices', ...
                'default', 0);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            number = obj.getParameter('number');
            if ~isempty(number) && number < 1
                valid = false;
                warning('Stack:InvalidNumber', ...
                    'Number of items must be at least 1');
            end
            
            % Check block size relative to total number
            blockSize = obj.getParameter('block-size');
            if ~isempty(blockSize) && ~isempty(number)
                if blockSize > number
                    warning('Stack:LargeBlockSize', ...
                        'Block size (%d) exceeds total items (%d)', ...
                        blockSize, number);
                end
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate stack operations
            
            fprintf('\n=== Stack Task Demo ===\n\n');
            
            % Create synthetic volume
            volume = zeros(256, 256, 50);
            for i = 1:50
                % Create varying circular patterns
                [X, Y] = meshgrid(-128:127, -128:127);
                R = sqrt(X.^2 + Y.^2);
                volume(:, :, i) = exp(-(R - 20*sin(i/10)).^2 / 500);
            end
            
            % Create task
            task = mufo.tasks.Stack();
            
            % Example 1: Basic stacking
            fprintf('1. Basic stacking operation:\n');
            task.setNumber(50);
            task.getStackInfo([256, 256], 50);
            
            % Example 2: Block processing for large data
            fprintf('\n2. Block processing for memory efficiency:\n');
            task.setNumber(1000);
            task.setBlockSize(100);
            task.getStackInfo([2048, 2048], 1000);
            
            % Example 3: Different stacking dimensions
            fprintf('\n3. Stacking along different dimensions:\n');
            task.setDimension(1);
            task.getStackInfo([256, 256], 50);
            
            % Example 4: Create stacking chain
            fprintf('\n4. Creating stacking chain:\n');
            chain = task.createStackingChain('slice_*.tif', 'volume.h5');
            fprintf('Chain created for stacking TIFF slices to HDF5\n');
            
            % Example 5: Visual demonstration
            fprintf('\n5. Visual demonstration:\n');
            task.demonstrateStacking(volume);
            
            fprintf('\nDemo complete.\n');
        end
    end
end