classdef DataLoader < handle
    % DataLoader - Unified data loading interface for MUFO
    %
    % This class provides a unified interface for loading various data formats
    % commonly used in tomographic reconstruction, including HDF5, TIFF, and raw files.
    %
    % Example:
    %   loader = mufo.data.DataLoader();
    %   data = loader.load('projections.h5', '/exchange/data');
    %   metadata = loader.getMetadata();
    
    properties (Access = private)
        currentFile         % Current file being loaded
        currentDataset      % Current dataset path (for HDF5)
        metadata            % Metadata structure
        dataInfo            % Information about loaded data
        cache               % Data cache for performance
        useCache            % Enable/disable caching
    end
    
    properties (Constant)
        SUPPORTED_FORMATS = {'.h5', '.hdf5', '.tif', '.tiff', '.raw', '.edf', '.fits'};
        HDF5_EXCHANGE_PATHS = struct(...
            'data', '/exchange/data', ...
            'dark', '/exchange/data_dark', ...
            'flat', '/exchange/data_white', ...
            'theta', '/exchange/theta');
    end
    
    events
        DataLoaded
        MetadataUpdated
    end
    
    methods
        function obj = DataLoader(varargin)
            % Constructor
            %
            % Parameters:
            %   'UseCache' - Enable caching (default: true)
            %   'CacheSize' - Maximum cache size in MB (default: 1000)
            
            p = inputParser;
            addParameter(p, 'UseCache', true, @islogical);
            addParameter(p, 'CacheSize', 1000, @isnumeric);
            parse(p, varargin{:});
            
            obj.useCache = p.Results.UseCache;
            obj.cache = containers.Map();
            obj.metadata = struct();
        end
        
        function data = load(obj, filename, varargin)
            % Load data from file
            %
            % Inputs:
            %   filename - Path to data file
            %   dataset - Dataset path for HDF5 (optional)
            %   'Range' - Range of slices to load [start, end]
            %   'Stride' - Step between slices
            %   'ROI' - Region of interest [x, y, width, height]
            %   'DataType' - Convert to specific data type
            %
            % Output:
            %   data - Loaded data array
            
            p = inputParser;
            addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
            addOptional(p, 'dataset', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Range', [], @isnumeric);
            addParameter(p, 'Stride', 1, @isnumeric);
            addParameter(p, 'ROI', [], @isnumeric);
            addParameter(p, 'DataType', '', @ischar);
            parse(p, filename, varargin{:});
            
            % Check if file exists
            if ~exist(filename, 'file')
                error('DataLoader:FileNotFound', 'File not found: %s', filename);
            end
            
            obj.currentFile = filename;
            obj.currentDataset = p.Results.dataset;
            
            % Check cache
            cacheKey = obj.getCacheKey(filename, p.Results);
            if obj.useCache && obj.cache.isKey(cacheKey)
                data = obj.cache(cacheKey);
                notify(obj, 'DataLoaded');
                return;
            end
            
            % Load based on file format
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case {'.h5', '.hdf5'}
                    data = obj.loadHDF5(filename, p.Results);
                    
                case {'.tif', '.tiff'}
                    data = obj.loadTIFF(filename, p.Results);
                    
                case '.raw'
                    data = obj.loadRAW(filename, p.Results);
                    
                case '.edf'
                    data = obj.loadEDF(filename, p.Results);
                    
                case '.fits'
                    data = obj.loadFITS(filename, p.Results);
                    
                otherwise
                    error('DataLoader:UnsupportedFormat', ...
                        'Unsupported file format: %s', ext);
            end
            
            % Apply data type conversion
            if ~isempty(p.Results.DataType)
                data = cast(data, p.Results.DataType);
            end
            
            % Update cache
            if obj.useCache
                obj.addToCache(cacheKey, data);
            end
            
            % Store data info
            obj.updateDataInfo(data);
            
            % Notify listeners
            notify(obj, 'DataLoaded');
        end
        
        function [projections, flats, darks, angles] = loadExchangeFormat(obj, filename)
            % Load complete dataset in HDF5 exchange format
            %
            % Input:
            %   filename - HDF5 file in exchange format
            %
            % Outputs:
            %   projections - Projection data
            %   flats - Flat field data
            %   darks - Dark field data
            %   angles - Projection angles
            
            if ~obj.isHDF5(filename)
                error('DataLoader:InvalidFormat', ...
                    'Exchange format requires HDF5 file');
            end
            
            % Load each component
            projections = obj.load(filename, obj.HDF5_EXCHANGE_PATHS.data);
            
            % Check for optional data
            info = h5info(filename);
            datasets = obj.listDatasets(info);
            
            flats = [];
            if any(strcmp(datasets, obj.HDF5_EXCHANGE_PATHS.flat))
                flats = obj.load(filename, obj.HDF5_EXCHANGE_PATHS.flat);
            end
            
            darks = [];
            if any(strcmp(datasets, obj.HDF5_EXCHANGE_PATHS.dark))
                darks = obj.load(filename, obj.HDF5_EXCHANGE_PATHS.dark);
            end
            
            angles = [];
            if any(strcmp(datasets, obj.HDF5_EXCHANGE_PATHS.theta))
                angles = obj.load(filename, obj.HDF5_EXCHANGE_PATHS.theta);
            end
            
            % Extract metadata
            obj.extractExchangeMetadata(filename);
        end
        
        function metadata = getMetadata(obj)
            % Get current metadata
            metadata = obj.metadata;
        end
        
        function info = getDataInfo(obj)
            % Get information about loaded data
            info = obj.dataInfo;
        end
        
        function previewData(obj, filename, varargin)
            % Preview data without fully loading
            %
            % Inputs:
            %   filename - Data file
            %   'Slice' - Slice number to preview
            %   'Downsample' - Downsampling factor
            
            p = inputParser;
            addRequired(p, 'filename');
            addParameter(p, 'Slice', 'middle', @(x) isnumeric(x) || ischar(x));
            addParameter(p, 'Downsample', 4, @isnumeric);
            parse(p, filename, varargin{:});
            
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case {'.h5', '.hdf5'}
                    obj.previewHDF5(filename, p.Results);
                    
                case {'.tif', '.tiff'}
                    obj.previewTIFF(filename, p.Results);
                    
                otherwise
                    % Load small portion
                    data = obj.load(filename, 'Range', [1, 1]);
                    obj.displayPreview(data);
            end
        end
        
        function convertFormat(obj, inputFile, outputFile, varargin)
            % Convert between data formats
            %
            % Inputs:
            %   inputFile - Input filename
            %   outputFile - Output filename
            %   'Compression' - Enable compression
            %   'ChunkSize' - Chunk size for processing
            
            p = inputParser;
            addRequired(p, 'inputFile');
            addRequired(p, 'outputFile');
            addParameter(p, 'Compression', true, @islogical);
            addParameter(p, 'ChunkSize', 100, @isnumeric);
            parse(p, inputFile, outputFile, varargin{:});
            
            % Determine output format
            [~, ~, outExt] = fileparts(outputFile);
            
            % Get input data info
            info = obj.getFileInfo(inputFile);
            
            % Process in chunks for large data
            nChunks = ceil(info.slices / p.Results.ChunkSize);
            
            writer = mufo.data.DataWriter('Compression', p.Results.Compression);
            
            h = waitbar(0, 'Converting format...');
            
            for i = 1:nChunks
                startIdx = (i-1) * p.Results.ChunkSize + 1;
                endIdx = min(i * p.Results.ChunkSize, info.slices);
                
                % Load chunk
                chunk = obj.load(inputFile, 'Range', [startIdx, endIdx]);
                
                % Write chunk
                if i == 1
                    writer.create(outputFile, size(chunk), class(chunk));
                end
                
                writer.writeSlices(chunk, startIdx);
                
                waitbar(i/nChunks, h);
            end
            
            close(h);
            
            % Copy metadata if applicable
            if strcmpi(outExt, '.h5') || strcmpi(outExt, '.hdf5')
                obj.copyMetadata(inputFile, outputFile);
            end
        end
    end
    
    methods (Access = private)
        function data = loadHDF5(obj, filename, options)
            % Load HDF5 data
            
            % Determine dataset path
            if isempty(options.dataset)
                % Try to find main dataset
                dataset = obj.findMainDataset(filename);
            else
                dataset = options.dataset;
            end
            
            if isempty(dataset)
                error('DataLoader:NoDataset', ...
                    'No dataset specified and none found automatically');
            end
            
            % Get dataset info
            info = h5info(filename, dataset);
            dims = info.Dataspace.Size;
            
            % Determine read parameters
            if isempty(options.Range)
                start = ones(1, length(dims));
                count = dims;
            else
                % Handle range specification
                start = ones(1, length(dims));
                count = dims;
                
                % Apply range to last dimension (usually slices)
                start(end) = options.Range(1);
                if length(options.Range) > 1
                    count(end) = options.Range(2) - options.Range(1) + 1;
                else
                    count(end) = 1;
                end
            end
            
            % Apply stride
            stride = ones(1, length(dims));
            if options.Stride > 1
                stride(end) = options.Stride;
                count(end) = ceil(count(end) / options.Stride);
            end
            
            % Read data
            data = h5read(filename, dataset, start, count, stride);
            
            % Apply ROI if specified
            if ~isempty(options.ROI) && length(options.ROI) == 4
                x = options.ROI(1);
                y = options.ROI(2);
                w = options.ROI(3);
                h = options.ROI(4);
                
                if ndims(data) == 3
                    data = data(y:y+h-1, x:x+w-1, :);
                else
                    data = data(y:y+h-1, x:x+w-1);
                end
            end
            
            % Extract metadata
            obj.extractHDF5Metadata(filename, dataset);
        end
        
        function data = loadTIFF(obj, filename, options)
            % Load TIFF data
            
            % Check if single file or pattern
            if contains(filename, '*') || contains(filename, '?')
                % File pattern
                files = dir(filename);
                if isempty(files)
                    error('DataLoader:NoFiles', ...
                        'No files match pattern: %s', filename);
                end
                
                % Sort files naturally
                [~, idx] = sort_nat({files.name});
                files = files(idx);
                
                % Determine range
                if isempty(options.Range)
                    fileRange = 1:length(files);
                else
                    fileRange = options.Range(1):options.Stride:min(options.Range(2), length(files));
                end
                
                % Read first image to get dimensions
                firstImg = imread(fullfile(files(1).folder, files(1).name));
                
                % Apply ROI to dimensions
                if ~isempty(options.ROI)
                    h = options.ROI(4);
                    w = options.ROI(3);
                else
                    [h, w] = size(firstImg);
                end
                
                % Preallocate
                data = zeros(h, w, length(fileRange), class(firstImg));
                
                % Load images
                for i = 1:length(fileRange)
                    img = imread(fullfile(files(fileRange(i)).folder, ...
                                         files(fileRange(i)).name));
                    
                    % Apply ROI
                    if ~isempty(options.ROI)
                        img = img(options.ROI(2):options.ROI(2)+options.ROI(4)-1, ...
                                 options.ROI(1):options.ROI(1)+options.ROI(3)-1);
                    end
                    
                    data(:, :, i) = img;
                end
                
            else
                % Single multi-page TIFF
                info = imfinfo(filename);
                nFrames = length(info);
                
                % Determine range
                if isempty(options.Range)
                    frameRange = 1:options.Stride:nFrames;
                else
                    frameRange = options.Range(1):options.Stride:min(options.Range(2), nFrames);
                end
                
                % Read first frame to get dimensions
                firstImg = imread(filename, 1);
                
                % Apply ROI to dimensions
                if ~isempty(options.ROI)
                    h = options.ROI(4);
                    w = options.ROI(3);
                else
                    [h, w] = size(firstImg);
                end
                
                % Preallocate
                data = zeros(h, w, length(frameRange), class(firstImg));
                
                % Load frames
                for i = 1:length(frameRange)
                    img = imread(filename, frameRange(i));
                    
                    % Apply ROI
                    if ~isempty(options.ROI)
                        img = img(options.ROI(2):options.ROI(2)+options.ROI(4)-1, ...
                                 options.ROI(1):options.ROI(1)+options.ROI(3)-1);
                    end
                    
                    data(:, :, i) = img;
                end
            end
            
            % Extract metadata
            obj.extractTIFFMetadata(filename);
        end
        
        function data = loadRAW(obj, filename, options)
            % Load raw binary data
            
            % Need metadata for raw files
            if ~isfield(obj.metadata, 'raw') || isempty(obj.metadata.raw)
                % Try to load metadata from accompanying file
                metaFile = [filename, '.info'];
                if exist(metaFile, 'file')
                    obj.metadata.raw = obj.loadRawMetadata(metaFile);
                else
                    error('DataLoader:NoMetadata', ...
                        'Raw file requires metadata (.info file or manual specification)');
                end
            end
            
            meta = obj.metadata.raw;
            
            % Open file
            fid = fopen(filename, 'rb');
            if fid == -1
                error('DataLoader:FileError', 'Cannot open file: %s', filename);
            end
            
            % Calculate dimensions and seek position
            if isempty(options.Range)
                sliceStart = 1;
                sliceCount = meta.slices;
            else
                sliceStart = options.Range(1);
                sliceCount = options.Range(2) - options.Range(1) + 1;
            end
            
            % Seek to start position
            bytesPerPixel = obj.getBytesPerPixel(meta.dataType);
            sliceBytes = meta.width * meta.height * bytesPerPixel;
            offset = (sliceStart - 1) * sliceBytes;
            fseek(fid, offset, 'bof');
            
            % Read data
            if options.Stride == 1
                % Read continuous block
                count = meta.width * meta.height * sliceCount;
                data = fread(fid, count, ['*', meta.dataType]);
                data = reshape(data, [meta.height, meta.width, sliceCount]);
            else
                % Read with stride
                nSlices = ceil(sliceCount / options.Stride);
                data = zeros(meta.height, meta.width, nSlices, meta.dataType);
                
                for i = 1:nSlices
                    slice = fread(fid, meta.width * meta.height, ['*', meta.dataType]);
                    data(:, :, i) = reshape(slice, [meta.height, meta.width]);
                    
                    % Skip to next slice
                    if i < nSlices
                        fseek(fid, (options.Stride - 1) * sliceBytes, 'cof');
                    end
                end
            end
            
            fclose(fid);
            
            % Apply ROI
            if ~isempty(options.ROI)
                data = data(options.ROI(2):options.ROI(2)+options.ROI(4)-1, ...
                           options.ROI(1):options.ROI(1)+options.ROI(3)-1, :);
            end
        end
        
        function data = loadEDF(obj, filename, options)
            % Load ESRF Data Format files
            
            % Read header
            header = obj.readEDFHeader(filename);
            
            % Open file
            fid = fopen(filename, 'rb');
            if fid == -1
                error('DataLoader:FileError', 'Cannot open file: %s', filename);
            end
            
            % Skip header
            fseek(fid, header.headerSize, 'bof');
            
            % Read data
            data = fread(fid, header.dim1 * header.dim2, ['*', header.dataType]);
            data = reshape(data, [header.dim2, header.dim1]);
            
            fclose(fid);
            
            % Handle byte order
            if ~strcmp(header.byteOrder, obj.getSystemByteOrder())
                data = swapbytes(data);
            end
            
            % Apply ROI
            if ~isempty(options.ROI)
                data = data(options.ROI(2):options.ROI(2)+options.ROI(4)-1, ...
                           options.ROI(1):options.ROI(1)+options.ROI(3)-1);
            end
            
            % Store metadata
            obj.metadata.edf = header;
        end
        
        function data = loadFITS(obj, filename, options)
            % Load FITS files
            
            % Use MATLAB's fitsread if available
            if exist('fitsread', 'file')
                data = fitsread(filename);
                info = fitsinfo(filename);
                obj.metadata.fits = info;
            else
                error('DataLoader:NoFITSSupport', ...
                    'FITS support requires Image Processing Toolbox');
            end
            
            % Apply ROI
            if ~isempty(options.ROI)
                if ndims(data) == 3
                    data = data(options.ROI(2):options.ROI(2)+options.ROI(4)-1, ...
                               options.ROI(1):options.ROI(1)+options.ROI(3)-1, :);
                else
                    data = data(options.ROI(2):options.ROI(2)+options.ROI(4)-1, ...
                               options.ROI(1):options.ROI(1)+options.ROI(3)-1);
                end
            end
            
            % Apply range for 3D data
            if ndims(data) == 3 && ~isempty(options.Range)
                data = data(:, :, options.Range(1):options.Stride:options.Range(2));
            end
        end
        
        function dataset = findMainDataset(obj, filename)
            % Find main dataset in HDF5 file
            
            info = h5info(filename);
            datasets = obj.listDatasets(info);
            
            % Look for common dataset names
            commonNames = {'/exchange/data', '/entry/data/data', ...
                          '/measurement/data', '/data'};
            
            for i = 1:length(commonNames)
                if any(strcmp(datasets, commonNames{i}))
                    dataset = commonNames{i};
                    return;
                end
            end
            
            % Return largest dataset
            if ~isempty(datasets)
                maxSize = 0;
                dataset = '';
                
                for i = 1:length(datasets)
                    dsInfo = h5info(filename, datasets{i});
                    dataSize = prod(dsInfo.Dataspace.Size);
                    if dataSize > maxSize
                        maxSize = dataSize;
                        dataset = datasets{i};
                    end
                end
            else
                dataset = '';
            end
        end
        
        function datasets = listDatasets(obj, info, prefix)
            % Recursively list all datasets in HDF5 structure
            
            if nargin < 3
                prefix = '';
            end
            
            datasets = {};
            
            % Add datasets at current level
            for i = 1:length(info.Datasets)
                datasets{end+1} = fullfile(prefix, info.Datasets(i).Name);
            end
            
            % Recursively check groups
            for i = 1:length(info.Groups)
                groupDatasets = obj.listDatasets(info.Groups(i), info.Groups(i).Name);
                datasets = [datasets, groupDatasets];
            end
        end
        
        function extractHDF5Metadata(obj, filename, dataset)
            % Extract metadata from HDF5 file
            
            obj.metadata.hdf5 = struct();
            
            % File-level attributes
            fileInfo = h5info(filename);
            obj.metadata.hdf5.file_attributes = obj.getAttributes(fileInfo);
            
            % Dataset attributes
            if ~isempty(dataset)
                dsInfo = h5info(filename, dataset);
                obj.metadata.hdf5.dataset_attributes = obj.getAttributes(dsInfo);
                obj.metadata.hdf5.datatype = dsInfo.Datatype;
                obj.metadata.hdf5.dimensions = dsInfo.Dataspace.Size;
            end
            
            % Look for specific metadata groups
            if any(strcmp(obj.listDatasets(fileInfo), '/measurement/instrument'))
                obj.metadata.instrument = obj.readHDF5Group(filename, '/measurement/instrument');
            end
            
            if any(strcmp(obj.listDatasets(fileInfo), '/measurement/sample'))
                obj.metadata.sample = obj.readHDF5Group(filename, '/measurement/sample');
            end
            
            notify(obj, 'MetadataUpdated');
        end
        
        function extractTIFFMetadata(obj, filename)
            % Extract metadata from TIFF file
            
            info = imfinfo(filename);
            
            obj.metadata.tiff = struct();
            obj.metadata.tiff.width = info(1).Width;
            obj.metadata.tiff.height = info(1).Height;
            obj.metadata.tiff.bitDepth = info(1).BitDepth;
            obj.metadata.tiff.compression = info(1).Compression;
            
            if length(info) > 1
                obj.metadata.tiff.numFrames = length(info);
            end
            
            % Extract any custom tags
            if isfield(info(1), 'UnknownTags')
                obj.metadata.tiff.customTags = info(1).UnknownTags;
            end
            
            notify(obj, 'MetadataUpdated');
        end
        
        function extractExchangeMetadata(obj, filename)
            % Extract metadata from HDF5 exchange format
            
            info = h5info(filename);
            
            % Standard exchange metadata locations
            metaPaths = {
                '/measurement/instrument/source/energy'
                '/measurement/instrument/detector/pixel_size'
                '/measurement/instrument/detector/distance'
                '/measurement/sample/name'
                '/process/acquisition/rotation/num_angles'
                '/process/acquisition/rotation/angle_range'
            };
            
            obj.metadata.exchange = struct();
            
            for i = 1:length(metaPaths)
                try
                    value = h5read(filename, metaPaths{i});
                    fieldName = obj.pathToFieldName(metaPaths{i});
                    obj.metadata.exchange.(fieldName) = value;
                catch
                    % Path doesn't exist
                end
            end
            
            notify(obj, 'MetadataUpdated');
        end
        
        function attrs = getAttributes(obj, info)
            % Extract attributes from HDF5 info structure
            
            attrs = struct();
            
            if isfield(info, 'Attributes')
                for i = 1:length(info.Attributes)
                    name = matlab.lang.makeValidName(info.Attributes(i).Name);
                    attrs.(name) = info.Attributes(i).Value;
                end
            end
        end
        
        function group = readHDF5Group(obj, filename, groupPath)
            % Read all datasets in an HDF5 group
            
            groupInfo = h5info(filename, groupPath);
            group = struct();
            
            % Read datasets
            for i = 1:length(groupInfo.Datasets)
                name = groupInfo.Datasets(i).Name;
                path = fullfile(groupPath, name);
                group.(name) = h5read(filename, path);
            end
            
            % Read attributes
            group.attributes = obj.getAttributes(groupInfo);
        end
        
        function header = readEDFHeader(obj, filename)
            % Read EDF header
            
            fid = fopen(filename, 'rb');
            if fid == -1
                error('DataLoader:FileError', 'Cannot open file: %s', filename);
            end
            
            % Read header (usually first 1024 or 512 bytes)
            headerData = fread(fid, 1024, '*char')';
            fclose(fid);
            
            % Parse header
            header = struct();
            header.headerSize = 1024;  % Default
            
            % Extract key-value pairs
            lines = strsplit(headerData, ';');
            for i = 1:length(lines)
                if contains(lines{i}, '=')
                    parts = strsplit(lines{i}, '=');
                    key = strtrim(parts{1});
                    value = strtrim(parts{2});
                    
                    switch lower(key)
                        case 'dim_1'
                            header.dim1 = str2double(value);
                        case 'dim_2'
                            header.dim2 = str2double(value);
                        case 'datatype'
                            header.dataType = obj.edfToMatlabType(value);
                        case 'byteorder'
                            header.byteOrder = value;
                        case 'size'
                            header.headerSize = str2double(value);
                    end
                end
            end
        end
        
        function matlabType = edfToMatlabType(obj, edfType)
            % Convert EDF data type to MATLAB type
            
            switch lower(edfType)
                case 'signedbyte'
                    matlabType = 'int8';
                case 'unsignedbyte'
                    matlabType = 'uint8';
                case 'signedshort'
                    matlabType = 'int16';
                case 'unsignedshort'
                    matlabType = 'uint16';
                case 'signedinteger'
                    matlabType = 'int32';
                case 'unsignedinteger'
                    matlabType = 'uint32';
                case 'floatvalue'
                    matlabType = 'single';
                case 'doublevalue'
                    matlabType = 'double';
                otherwise
                    matlabType = 'single';
            end
        end
        
        function meta = loadRawMetadata(obj, metaFile)
            % Load metadata for raw files
            
            if endsWith(metaFile, '.json')
                jsonStr = fileread(metaFile);
                meta = jsondecode(jsonStr);
            else
                % Simple text format
                meta = struct();
                
                fid = fopen(metaFile, 'r');
                while ~feof(fid)
                    line = fgetl(fid);
                    if contains(line, '=')
                        parts = strsplit(line, '=');
                        key = strtrim(parts{1});
                        value = strtrim(parts{2});
                        
                        % Try to convert to number
                        numVal = str2double(value);
                        if ~isnan(numVal)
                            meta.(key) = numVal;
                        else
                            meta.(key) = value;
                        end
                    end
                end
                fclose(fid);
            end
        end
        
        function updateDataInfo(obj, data)
            % Update data information structure
            
            obj.dataInfo = struct();
            obj.dataInfo.size = size(data);
            obj.dataInfo.ndims = ndims(data);
            obj.dataInfo.class = class(data);
            obj.dataInfo.min = min(data(:));
            obj.dataInfo.max = max(data(:));
            obj.dataInfo.mean = mean(data(:), 'omitnan');
            obj.dataInfo.std = std(data(:), 'omitnan');
            obj.dataInfo.memoryMB = numel(data) * obj.getBytesPerPixel(class(data)) / 1e6;
            
            % Check for special values
            obj.dataInfo.hasNaN = any(isnan(data(:)));
            obj.dataInfo.hasInf = any(isinf(data(:)));
            obj.dataInfo.numNaN = sum(isnan(data(:)));
            obj.dataInfo.numInf = sum(isinf(data(:)));
        end
        
        function info = getFileInfo(obj, filename)
            % Get file information without loading data
            
            info = struct();
            info.filename = filename;
            
            [~, ~, ext] = fileparts(filename);
            info.format = ext;
            
            fileInfo = dir(filename);
            if ~isempty(fileInfo)
                info.bytes = fileInfo.bytes;
                info.date = fileInfo.date;
            end
            
            switch lower(ext)
                case {'.h5', '.hdf5'}
                    h5info_data = h5info(filename);
                    dataset = obj.findMainDataset(filename);
                    if ~isempty(dataset)
                        dsInfo = h5info(filename, dataset);
                        info.dimensions = dsInfo.Dataspace.Size;
                        info.datatype = dsInfo.Datatype.Class;
                        
                        if length(info.dimensions) >= 3
                            info.slices = info.dimensions(3);
                        else
                            info.slices = 1;
                        end
                    end
                    
                case {'.tif', '.tiff'}
                    if contains(filename, '*') || contains(filename, '?')
                        files = dir(filename);
                        info.numFiles = length(files);
                        info.slices = info.numFiles;
                    else
                        imInfo = imfinfo(filename);
                        info.width = imInfo(1).Width;
                        info.height = imInfo(1).Height;
                        info.slices = length(imInfo);
                        info.bitDepth = imInfo(1).BitDepth;
                    end
            end
        end
        
        function fieldName = pathToFieldName(obj, path)
            % Convert HDF5 path to valid field name
            
            parts = strsplit(path, '/');
            parts(cellfun(@isempty, parts)) = [];
            fieldName = strjoin(parts, '_');
            fieldName = matlab.lang.makeValidName(fieldName);
        end
        
        function tf = isHDF5(obj, filename)
            % Check if file is HDF5
            [~, ~, ext] = fileparts(filename);
            tf = any(strcmpi(ext, {'.h5', '.hdf5'}));
        end
        
        function byteOrder = getSystemByteOrder(obj)
            % Get system byte order
            [~, ~, endian] = computer;
            if endian == 'L'
                byteOrder = 'littleEndian';
            else
                byteOrder = 'bigEndian';
            end
        end
        
        function bytes = getBytesPerPixel(obj, dataType)
            % Get bytes per pixel for data type
            
            switch dataType
                case {'int8', 'uint8'}
                    bytes = 1;
                case {'int16', 'uint16'}
                    bytes = 2;
                case {'int32', 'uint32', 'single'}
                    bytes = 4;
                case {'int64', 'uint64', 'double'}
                    bytes = 8;
                otherwise
                    bytes = 4;
            end
        end
        
        function key = getCacheKey(obj, filename, options)
            % Generate cache key
            
            keyParts = {filename};
            
            if ~isempty(options.dataset)
                keyParts{end+1} = options.dataset;
            end
            
            if ~isempty(options.Range)
                keyParts{end+1} = sprintf('range_%d_%d', options.Range);
            end
            
            if options.Stride > 1
                keyParts{end+1} = sprintf('stride_%d', options.Stride);
            end
            
            if ~isempty(options.ROI)
                keyParts{end+1} = sprintf('roi_%d_%d_%d_%d', options.ROI);
            end
            
            key = strjoin(keyParts, '_');
        end
        
        function addToCache(obj, key, data)
            % Add data to cache with size management
            
            % Simple size-based cache management
            maxCacheSize = 1000 * 1e6;  % 1GB in bytes
            dataSize = numel(data) * obj.getBytesPerPixel(class(data));
            
            % Clear cache if needed
            if obj.cache.Count > 0
                currentSize = 0;
                keys = obj.cache.keys;
                for i = 1:length(keys)
                    cachedData = obj.cache(keys{i});
                    currentSize = currentSize + numel(cachedData) * ...
                                 obj.getBytesPerPixel(class(cachedData));
                end
                
                if currentSize + dataSize > maxCacheSize
                    % Clear oldest entries (simple FIFO)
                    remove(obj.cache, keys{1});
                end
            end
            
            obj.cache(key) = data;
        end
        
        function displayPreview(obj, data)
            % Display data preview
            
            figure('Name', 'Data Preview');
            
            if ndims(data) == 3
                % Show middle slice
                midSlice = round(size(data, 3) / 2);
                imagesc(data(:, :, midSlice));
                title(sprintf('Slice %d / %d', midSlice, size(data, 3)));
            else
                imagesc(data);
                title('2D Data');
            end
            
            colormap gray;
            axis image;
            colorbar;
            
            % Add info
            xlabel(sprintf('Size: %s | Type: %s | Range: [%.2f, %.2f]', ...
                mat2str(size(data)), class(data), min(data(:)), max(data(:))));
        end
        
        function copyMetadata(obj, sourceFile, destFile)
            % Copy metadata between files
            
            % For HDF5 files
            if obj.isHDF5(sourceFile) && obj.isHDF5(destFile)
                % Copy attributes
                sourceInfo = h5info(sourceFile);
                
                % Copy file-level attributes
                for i = 1:length(sourceInfo.Attributes)
                    attr = sourceInfo.Attributes(i);
                    try
                        h5writeatt(destFile, '/', attr.Name, attr.Value);
                    catch
                        % Attribute might already exist or be write-protected
                    end
                end
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate DataLoader capabilities
            
            fprintf('\n=== DataLoader Demo ===\n\n');
            
            % Create loader
            loader = mufo.data.DataLoader();
            
            % Example 1: Load HDF5 exchange format
            fprintf('1. Loading HDF5 exchange format:\n');
            try
                filename = 'example_data.h5';
                if exist(filename, 'file')
                    [proj, flat, dark, angles] = loader.loadExchangeFormat(filename);
                    fprintf('   Loaded: projections %s, flats %s, darks %s\n', ...
                        mat2str(size(proj)), mat2str(size(flat)), mat2str(size(dark)));
                else
                    fprintf('   Example file not found\n');
                end
            catch ME
                fprintf('   Error: %s\n', ME.message);
            end
            
            % Example 2: Load with options
            fprintf('\n2. Loading with ROI and range:\n');
            try
                data = loader.load('test.tif', ...
                    'Range', [1, 10], ...
                    'ROI', [100, 100, 512, 512]);
                fprintf('   Loaded data: %s\n', mat2str(size(data)));
            catch ME
                fprintf('   Error: %s\n', ME.message);
            end
            
            % Example 3: Format conversion
            fprintf('\n3. Format conversion example:\n');
            fprintf('   Would convert TIFF stack to HDF5...\n');
            
            fprintf('\nDemo complete.\n');
        end
    end
end