classdef Write < mufo.core.UFOTask
    % Write - UFO write task wrapper with format support
    %
    % This task writes data to various formats including HDF5, TIFF,
    % and raw files. It supports compression and data type conversion.
    %
    % Example:
    %   task = mufo.tasks.Write();
    %   task.setFilename('output/slice-%04d.tif');
    %   task.setBits(16);  % Convert to 16-bit
    
    properties (Access = private)
        outputInfo      % Information about output configuration
    end
    
    methods
        function obj = Write(varargin)
            % Constructor
            obj@mufo.core.UFOTask('write');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setFilename(obj, filename)
            % Set output filename with format validation
            %
            % Input:
            %   filename - Output filename or pattern
            
            obj.setParameter('filename', filename);
            
            % Auto-detect format
            [path, name, ext] = fileparts(filename);
            
            % Create output directory if needed
            if ~isempty(path) && ~exist(path, 'dir')
                mkdir(path);
            end
            
            % Validate format
            supportedFormats = {'.tif', '.tiff', '.h5', '.hdf5', '.raw', '.edf'};
            if ~any(strcmpi(ext, supportedFormats))
                warning('Write:UnsupportedFormat', ...
                    'Format %s may not be supported by UFO', ext);
            end
            
            % Store output info
            obj.outputInfo.path = path;
            obj.outputInfo.pattern = name;
            obj.outputInfo.format = ext;
        end
        
        function setBits(obj, bits)
            % Set output bit depth
            %
            % Input:
            %   bits - Bit depth (8, 16, or 32)
            
            if ~any(bits == [8, 16, 32])
                error('Write:InvalidBits', ...
                    'Bits must be 8, 16, or 32');
            end
            
            obj.setParameter('bits', bits);
        end
        
        function setCompression(obj, enabled, level)
            % Set compression options
            %
            % Inputs:
            %   enabled - Enable compression (logical)
            %   level - Compression level (1-9, optional)
            
            if enabled
                obj.setParameter('compression', true);
                
                if nargin > 2 && ~isempty(level)
                    if level < 1 || level > 9
                        error('Write:InvalidCompression', ...
                            'Compression level must be between 1 and 9');
                    end
                    obj.setParameter('compression-level', level);
                end
            else
                obj.setParameter('compression', false);
            end
        end
        
        function setHDF5Options(obj, dataset, chunks, compression)
            % Set HDF5-specific options
            %
            % Inputs:
            %   dataset - Dataset name (default: '/exchange/data')
            %   chunks - Chunk size [x, y, z] (optional)
            %   compression - Compression level (optional)
            
            if nargin < 2 || isempty(dataset)
                dataset = '/exchange/data';
            end
            
            obj.setParameter('hdf5-dataset', dataset);
            
            if nargin > 2 && ~isempty(chunks)
                % UFO expects comma-separated string
                chunkStr = sprintf('%d,%d,%d', chunks(1), chunks(2), chunks(3));
                obj.setParameter('hdf5-chunk-size', chunkStr);
            end
            
            if nargin > 3 && ~isempty(compression)
                obj.setCompression(true, compression);
            end
        end
        
        function setROI(obj, x, y, width, height)
            % Set region of interest to write
            %
            % Inputs:
            %   x, y - Starting position
            %   width, height - Region size
            
            obj.setParameter('x', x);
            obj.setParameter('y', y);
            obj.setParameter('width', width);
            obj.setParameter('height', height);
        end
        
        function setRescale(obj, minVal, maxVal)
            % Set rescaling range
            %
            % Inputs:
            %   minVal - Minimum value for rescaling
            %   maxVal - Maximum value for rescaling
            
            obj.setParameter('minimum', minVal);
            obj.setParameter('maximum', maxVal);
            obj.setParameter('rescale', true);
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get output settings
                outputDir = config.get('output.directory', 'output');
                format = config.get('output.format', 'tiff');
                dtype = config.get('output.dtype', 'float32');
                
                % Build filename pattern
                if strcmpi(format, 'tiff') || strcmpi(format, 'tif')
                    filename = fullfile(outputDir, 'slice_%05d.tif');
                elseif strcmpi(format, 'hdf5') || strcmpi(format, 'h5')
                    filename = fullfile(outputDir, 'reconstructed.h5');
                else
                    filename = fullfile(outputDir, 'output');
                end
                
                obj.setFilename(filename);
                
                % Set bit depth based on dtype
                switch lower(dtype)
                    case {'uint8', 'int8'}
                        obj.setBits(8);
                    case {'uint16', 'int16'}
                        obj.setBits(16);
                    case {'float32', 'single', 'uint32', 'int32'}
                        obj.setBits(32);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'filename')
                    obj.setFilename(config.filename);
                end
                
                if isfield(config, 'bits')
                    obj.setBits(config.bits);
                end
                
                if isfield(config, 'compression')
                    level = [];
                    if isfield(config, 'compressionLevel')
                        level = config.compressionLevel;
                    end
                    obj.setCompression(config.compression, level);
                end
                
                if isfield(config, 'rescale') && config.rescale
                    obj.setRescale(config.minimum, config.maximum);
                end
            end
        end
        
        function validateOutput(obj)
            % Validate output configuration
            filename = obj.getParameter('filename');
            
            if isempty(filename)
                error('Write:NoFilename', 'No output filename specified');
            end
            
            % Check if output directory is writable
            [path, ~, ~] = fileparts(filename);
            if isempty(path)
                path = '.';
            end
            
            if ~exist(path, 'dir')
                try
                    mkdir(path);
                catch ME
                    error('Write:DirectoryError', ...
                        'Cannot create output directory: %s', ME.message);
                end
            end
            
            % Test write permissions
            testFile = fullfile(path, '.write_test');
            try
                fid = fopen(testFile, 'w');
                if fid == -1
                    error('Write:NoPermission', ...
                        'No write permission in output directory');
                end
                fclose(fid);
                delete(testFile);
            catch ME
                error('Write:WriteError', ...
                    'Cannot write to output directory: %s', ME.message);
            end
        end
        
        function estimateSize(obj, dims, dtype)
            % Estimate output file size
            %
            % Inputs:
            %   dims - Data dimensions [x, y, z]
            %   dtype - Data type (optional)
            % Output:
            %   sizeInfo - Structure with size information
            
            if nargin < 3
                bits = obj.getParameter('bits');
                if isempty(bits)
                    bits = 32;  % Default
                end
                bytesPerPixel = bits / 8;
            else
                switch dtype
                    case {'uint8', 'int8'}
                        bytesPerPixel = 1;
                    case {'uint16', 'int16'}
                        bytesPerPixel = 2;
                    case {'float32', 'single', 'uint32', 'int32'}
                        bytesPerPixel = 4;
                    case {'double', 'float64'}
                        bytesPerPixel = 8;
                    otherwise
                        bytesPerPixel = 4;
                end
            end
            
            % Calculate uncompressed size
            totalPixels = prod(dims);
            uncompressedSize = totalPixels * bytesPerPixel;
            
            % Estimate compressed size
            compression = obj.getParameter('compression');
            if compression
                % Rough estimates based on typical compression ratios
                format = obj.outputInfo.format;
                switch lower(format)
                    case {'.tif', '.tiff'}
                        compressionRatio = 0.7;  % LZW typically
                    case {'.h5', '.hdf5'}
                        compressionRatio = 0.5;  % Better compression
                    otherwise
                        compressionRatio = 1.0;  % No compression
                end
                compressedSize = uncompressedSize * compressionRatio;
            else
                compressedSize = uncompressedSize;
            end
            
            % Create info structure
            sizeInfo = struct();
            sizeInfo.dims = dims;
            sizeInfo.bytesPerPixel = bytesPerPixel;
            sizeInfo.uncompressedBytes = uncompressedSize;
            sizeInfo.compressedBytes = compressedSize;
            sizeInfo.uncompressedGB = uncompressedSize / 1e9;
            sizeInfo.compressedGB = compressedSize / 1e9;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\nEstimated output size:\n');
                fprintf('  Dimensions: %dx%dx%d\n', dims(1), dims(2), dims(3));
                fprintf('  Bytes per pixel: %d\n', bytesPerPixel);
                fprintf('  Uncompressed: %.2f GB\n', sizeInfo.uncompressedGB);
                if compression
                    fprintf('  Compressed: %.2f GB (estimated)\n', sizeInfo.compressedGB);
                end
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define write task parameters
            
            % Filename (required)
            obj.addParameterDef('filename', 'string', ...
                'description', 'Output filename or pattern', ...
                'required', true);
            
            % Data type
            obj.addParameterDef('bits', 'integer', ...
                'description', 'Output bit depth', ...
                'values', [8, 16, 32], ...
                'default', 32);
            
            % Compression
            obj.addParameterDef('compression', 'logical', ...
                'description', 'Enable compression', ...
                'default', false);
            
            obj.addParameterDef('compression-level', 'integer', ...
                'description', 'Compression level', ...
                'range', [1, 9], ...
                'default', 6);
            
            % Rescaling
            obj.addParameterDef('rescale', 'logical', ...
                'description', 'Enable rescaling', ...
                'default', false);
            
            obj.addParameterDef('minimum', 'numeric', ...
                'description', 'Minimum value for rescaling');
            
            obj.addParameterDef('maximum', 'numeric', ...
                'description', 'Maximum value for rescaling');
            
            % Region of interest
            obj.addParameterDef('x', 'integer', ...
                'description', 'X offset', ...
                'range', [0, inf]);
            
            obj.addParameterDef('y', 'integer', ...
                'description', 'Y offset', ...
                'range', [0, inf]);
            
            obj.addParameterDef('width', 'integer', ...
                'description', 'Width to write', ...
                'range', [1, inf]);
            
            obj.addParameterDef('height', 'integer', ...
                'description', 'Height to write', ...
                'range', [1, inf]);
            
            % HDF5 specific
            obj.addParameterDef('hdf5-dataset', 'string', ...
                'description', 'HDF5 dataset name', ...
                'default', '/exchange/data');
            
            obj.addParameterDef('hdf5-chunk-size', 'string', ...
                'description', 'HDF5 chunk dimensions');
            
            % TIFF specific
            obj.addParameterDef('tiff-compression', 'string', ...
                'description', 'TIFF compression method', ...
                'values', {'none', 'lzw', 'packbits'});
            
            % Overwrite
            obj.addParameterDef('overwrite', 'logical', ...
                'description', 'Overwrite existing files', ...
                'default', true);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            % Check if filename is set
            if isempty(obj.getParameter('filename'))
                valid = false;
                return;
            end
            
            % Validate rescaling parameters
            if obj.getParameter('rescale')
                minVal = obj.getParameter('minimum');
                maxVal = obj.getParameter('maximum');
                
                if isempty(minVal) || isempty(maxVal)
                    valid = false;
                    warning('Write:MissingRescale', ...
                        'Rescaling enabled but min/max values not set');
                elseif minVal >= maxVal
                    valid = false;
                    warning('Write:InvalidRescale', ...
                        'Minimum value must be less than maximum');
                end
            end
            
            % Validate region parameters
            x = obj.getParameter('x');
            y = obj.getParameter('y');
            width = obj.getParameter('width');
            height = obj.getParameter('height');
            
            hasRegion = ~isempty(x) || ~isempty(y) || ...
                       ~isempty(width) || ~isempty(height);
            
            if hasRegion
                % If any region parameter is set, all must be set
                if isempty(x) || isempty(y) || ...
                   isempty(width) || isempty(height)
                    valid = false;
                    warning('Write:IncompleteRegion', ...
                        'All region parameters must be set');
                end
            end
        end
    end
end