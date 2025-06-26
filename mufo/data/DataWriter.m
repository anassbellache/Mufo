classdef DataWriter < handle
    % DataWriter - Unified data writing interface for MUFO
    %
    % This class provides a unified interface for writing data to various formats
    % with support for compression, chunking, and metadata.
    %
    % Example:
    %   writer = mufo.data.DataWriter();
    %   writer.write(data, 'output.h5', '/exchange/data');
    
    properties (Access = private)
        currentFile         % Current output file
        fileHandle         % File handle for incremental writing
        compression        % Compression settings
        chunkSize          % Chunk size for HDF5
        metadata           % Metadata to write
    end
    
    properties (Constant)
        SUPPORTED_FORMATS = {'.h5', '.hdf5', '.tif', '.tiff', '.raw'};
        DEFAULT_COMPRESSION = struct(...
            'h5', 'gzip', ...
            'h5_level', 4, ...
            'tiff', 'lzw');
    end
    
    events
        DataWritten
        FileCreated
        FileClosed
    end
    
    methods
        function obj = DataWriter(varargin)
            % Constructor
            %
            % Parameters:
            %   'Compression' - Enable compression (default: true)
            %   'ChunkSize' - Chunk size for HDF5 [x y z]
            
            p = inputParser;
            addParameter(p, 'Compression', true, @islogical);
            addParameter(p, 'ChunkSize', [64, 64, 64], @isnumeric);
            parse(p, varargin{:});
            
            obj.compression = p.Results.Compression;
            obj.chunkSize = p.Results.ChunkSize;
            obj.metadata = struct();
        end
        
        function write(obj, data, filename, varargin)
            % Write data to file
            %
            % Inputs:
            %   data - Data array to write
            %   filename - Output filename
            %   dataset - HDF5 dataset path (optional)
            %   'Overwrite' - Overwrite existing file
            %   'Append' - Append to existing file
            %   'Metadata' - Metadata structure to include
            
            p = inputParser;
            addRequired(p, 'data', @isnumeric);
            addRequired(p, 'filename', @ischar);
            addOptional(p, 'dataset', '/data', @ischar);
            addParameter(p, 'Overwrite', true, @islogical);
            addParameter(p, 'Append', false, @islogical);
            addParameter(p, 'Metadata', struct(), @isstruct);
            parse(p, data, filename, varargin{:});
            
            % Check file format
            [~, ~, ext] = fileparts(filename);
            if ~any(strcmpi(ext, obj.SUPPORTED_FORMATS))
                error('DataWriter:UnsupportedFormat', ...
                    'Unsupported format: %s', ext);
            end
            
            % Check if file exists
            if exist(filename, 'file') && ~p.Results.Overwrite && ~p.Results.Append
                error('DataWriter:FileExists', ...
                    'File already exists: %s', filename);
            end
            
            % Merge metadata
            obj.metadata = obj.mergeStructs(obj.metadata, p.Results.Metadata);
            
            % Write based on format
            switch lower(ext)
                case {'.h5', '.hdf5'}
                    obj.writeHDF5(data, filename, p.Results.dataset, p.Results);
                    
                case {'.tif', '.tiff'}
                    obj.writeTIFF(data, filename, p.Results);
                    
                case '.raw'
                    obj.writeRAW(data, filename, p.Results);
                    
                otherwise
                    error('DataWriter:UnsupportedFormat', ...
                        'Format not implemented: %s', ext);
            end
            
            notify(obj, 'DataWritten');
        end
        
        function create(obj, filename, dimensions, datatype, varargin)
            % Create file for incremental writing
            %
            % Inputs:
            %   filename - Output filename
            %   dimensions - Data dimensions [x y z]
            %   datatype - Data type
            %   'Dataset' - HDF5 dataset path
            
            p = inputParser;
            addRequired(p, 'filename');
            addRequired(p, 'dimensions');
            addRequired(p, 'datatype');
            addParameter(p, 'Dataset', '/data', @ischar);
            parse(p, filename, dimensions, datatype, varargin{:});
            
            obj.currentFile = filename;
            
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case {'.h5', '.hdf5'}
                    obj.createHDF5(filename, dimensions, datatype, p.Results.Dataset);
                    
                case {'.tif', '.tiff'}
                    % TIFF doesn't need pre-creation
                    obj.fileHandle = [];
                    
                case '.raw'
                    obj.createRAW(filename, dimensions, datatype);
                    
                otherwise
                    error('DataWriter:UnsupportedFormat', ...
                        'Format not supported for incremental writing');
            end
            
            notify(obj, 'FileCreated');
        end
        
        function writeSlices(obj, data, startSlice)
            % Write slices incrementally
            %
            % Inputs:
            %   data - Slice data (2D or 3D)
            %   startSlice - Starting slice index
            
            if isempty(obj.currentFile)
                error('DataWriter:NoFile', ...
                    'No file created. Call create() first.');
            end
            
            [~, ~, ext] = fileparts(obj.currentFile);
            
            switch lower(ext)
                case {'.h5', '.hdf5'}
                    obj.writeHDF5Slices(data, startSlice);
                    
                case {'.tif', '.tiff'}
                    obj.writeTIFFSlices(data, startSlice);
                    
                case '.raw'
                    obj.writeRAWSlices(data, startSlice);
            end
        end
        
        function close(obj)
            % Close current file
            
            if ~isempty(obj.fileHandle)
                % Close file handle based on type
                if isa(obj.fileHandle, 'H5ML.id')
                    H5F.close(obj.fileHandle);
                elseif isnumeric(obj.fileHandle)
                    fclose(obj.fileHandle);
                end
                
                obj.fileHandle = [];
            end
            
            obj.currentFile = '';
            notify(obj, 'FileClosed');
        end
        
        function writeExchangeFormat(obj, filename, projections, varargin)
            % Write complete dataset in HDF5 exchange format
            %
            % Inputs:
            %   filename - Output HDF5 file
            %   projections - Projection data
            %   'Flats' - Flat field data
            %   'Darks' - Dark field data
            %   'Angles' - Projection angles
            %   'Metadata' - Additional metadata
            
            p = inputParser;
            addRequired(p, 'filename');
            addRequired(p, 'projections');
            addParameter(p, 'Flats', [], @isnumeric);
            addParameter(p, 'Darks', [], @isnumeric);
            addParameter(p, 'Angles', [], @isnumeric);
            addParameter(p, 'Metadata', struct(), @isstruct);
            parse(p, filename, projections, varargin{:});
            
            % Create file
            if exist(filename, 'file')
                delete(filename);
            end
            
            % Write projections
            obj.write(projections, filename, '/exchange/data');
            
            % Write flats if provided
            if ~isempty(p.Results.Flats)
                obj.write(p.Results.Flats, filename, '/exchange/data_white', ...
                    'Append', true);
            end
            
            % Write darks if provided
            if ~isempty(p.Results.Darks)
                obj.write(p.Results.Darks, filename, '/exchange/data_dark', ...
                    'Append', true);
            end
            
            % Write angles
            if ~isempty(p.Results.Angles)
                obj.write(p.Results.Angles, filename, '/exchange/theta', ...
                    'Append', true);
            else
                % Generate default angles
                nAngles = size(projections, 1);
                angles = linspace(0, 2*pi, nAngles+1);
                angles(end) = [];
                obj.write(angles(:), filename, '/exchange/theta', ...
                    'Append', true);
            end
            
            % Write metadata
            obj.writeExchangeMetadata(filename, p.Results.Metadata);
        end
        
        function setMetadata(obj, metadata)
            % Set metadata to be written
            obj.metadata = metadata;
        end
        
        function addMetadata(obj, key, value)
            % Add metadata entry
            obj.metadata.(key) = value;
        end
    end
    
    methods (Access = private)
        function writeHDF5(obj, data, filename, dataset, options)
            % Write HDF5 file
            
            % Determine write mode
            if options.Append && exist(filename, 'file')
                % Append mode - dataset must not exist
                try
                    h5info(filename, dataset);
                    error('DataWriter:DatasetExists', ...
                        'Dataset already exists: %s', dataset);
                catch ME
                    if ~contains(ME.identifier, 'DatasetNotFound')
                        rethrow(ME);
                    end
                end
            elseif exist(filename, 'file') && options.Overwrite
                % Try to delete dataset if it exists
                try
                    h5info(filename, dataset);
                    % Dataset exists, need to delete file and recreate
                    % (HDF5 doesn't support dataset deletion easily)
                    delete(filename);
                catch
                    % Dataset doesn't exist, OK to proceed
                end
            end
            
            % Create dataset with compression
            if obj.compression
                h5create(filename, dataset, size(data), ...
                    'Datatype', class(data), ...
                    'ChunkSize', obj.calculateChunkSize(size(data)), ...
                    'Deflate', obj.DEFAULT_COMPRESSION.h5_level);
            else
                h5create(filename, dataset, size(data), ...
                    'Datatype', class(data));
            end
            
            % Write data
            h5write(filename, dataset, data);
            
            % Write attributes
            obj.writeHDF5Attributes(filename, dataset);
        end
        
        function writeTIFF(obj, data, filename, options)
            % Write TIFF file(s)
            
            if ndims(data) == 3
                % Multi-page TIFF or image sequence
                [h, w, n] = size(data);
                
                if contains(filename, '%')
                    % Image sequence
                    [path, name, ext] = fileparts(filename);
                    
                    for i = 1:n
                        fname = fullfile(path, sprintf(name, i));
                        fname = [fname, ext];
                        
                        if obj.compression
                            imwrite(data(:,:,i), fname, ...
                                'Compression', obj.DEFAULT_COMPRESSION.tiff);
                        else
                            imwrite(data(:,:,i), fname, ...
                                'Compression', 'none');
                        end
                    end
                else
                    % Multi-page TIFF
                    for i = 1:n
                        if i == 1
                            mode = 'overwrite';
                        else
                            mode = 'append';
                        end
                        
                        if obj.compression
                            imwrite(data(:,:,i), filename, ...
                                'WriteMode', mode, ...
                                'Compression', obj.DEFAULT_COMPRESSION.tiff);
                        else
                            imwrite(data(:,:,i), filename, ...
                                'WriteMode', mode, ...
                                'Compression', 'none');
                        end
                    end
                end
            else
                % Single image
                if obj.compression
                    imwrite(data, filename, ...
                        'Compression', obj.DEFAULT_COMPRESSION.tiff);
                else
                    imwrite(data, filename, ...
                        'Compression', 'none');
                end
            end
        end
        
        function writeRAW(obj, data, filename, options)
            % Write raw binary file
            
            % Write data
            fid = fopen(filename, 'wb');
            if fid == -1
                error('DataWriter:FileError', ...
                    'Cannot create file: %s', filename);
            end
            
            fwrite(fid, data, class(data));
            fclose(fid);
            
            % Write metadata file
            metaFile = [filename, '.info'];
            obj.writeRAWMetadata(metaFile, data);
        end
        
        function createHDF5(obj, filename, dims, datatype, dataset)
            % Create HDF5 file for incremental writing
            
            % Create dataset
            if obj.compression
                h5create(filename, dataset, dims, ...
                    'Datatype', datatype, ...
                    'ChunkSize', obj.calculateChunkSize(dims), ...
                    'Deflate', obj.DEFAULT_COMPRESSION.h5_level);
            else
                h5create(filename, dataset, dims, ...
                    'Datatype', datatype);
            end
            
            % Store info for incremental writing
            obj.fileHandle = struct();
            obj.fileHandle.dataset = dataset;
            obj.fileHandle.dims = dims;
        end
        
        function createRAW(obj, filename, dims, datatype)
            % Create RAW file for incremental writing
            
            % Open file
            fid = fopen(filename, 'wb');
            if fid == -1
                error('DataWriter:FileError', ...
                    'Cannot create file: %s', filename);
            end
            
            obj.fileHandle = fid;
            
            % Write metadata
            metaFile = [filename, '.info'];
            meta = struct();
            meta.width = dims(1);
            meta.height = dims(2);
            if length(dims) > 2
                meta.slices = dims(3);
            else
                meta.slices = 1;
            end
            meta.dataType = datatype;
            obj.writeRAWMetadata(metaFile, [], meta);
        end
        
        function writeHDF5Slices(obj, data, startSlice)
            % Write slices to HDF5
            
            if ndims(data) == 2
                data = reshape(data, [size(data), 1]);
            end
            
            nSlices = size(data, 3);
            
            % Determine write position
            start = [1, 1, startSlice];
            count = [size(data, 1), size(data, 2), nSlices];
            
            % Write to dataset
            h5write(obj.currentFile, obj.fileHandle.dataset, data, start, count);
        end
        
        function writeTIFFSlices(obj, data, startSlice)
            % Write slices to TIFF
            
            if ndims(data) == 2
                data = reshape(data, [size(data), 1]);
            end
            
            nSlices = size(data, 3);
            
            % Append slices
            for i = 1:nSlices
                if obj.compression
                    imwrite(data(:,:,i), obj.currentFile, ...
                        'WriteMode', 'append', ...
                        'Compression', obj.DEFAULT_COMPRESSION.tiff);
                else
                    imwrite(data(:,:,i), obj.currentFile, ...
                        'WriteMode', 'append', ...
                        'Compression', 'none');
                end
            end
        end
        
        function writeRAWSlices(obj, data, startSlice)
            % Write slices to RAW file
            
            if obj.fileHandle == -1
                error('DataWriter:InvalidHandle', 'File not open');
            end
            
            % Write data
            fwrite(obj.fileHandle, data, class(data));
        end
        
        function chunks = calculateChunkSize(obj, dims)
            % Calculate optimal chunk size for HDF5
            
            % Use provided chunk size or calculate
            if length(obj.chunkSize) == length(dims)
                chunks = min(obj.chunkSize, dims);
            else
                % Auto-calculate
                chunks = dims;
                
                % For 3D data, chunk by slices
                if length(dims) == 3
                    chunks(3) = min(10, dims(3));
                end
                
                % Limit chunk size for memory efficiency
                maxChunkElements = 1e6;  % 1M elements per chunk
                while prod(chunks) > maxChunkElements
                    [~, maxDim] = max(chunks);
                    chunks(maxDim) = ceil(chunks(maxDim) / 2);
                end
            end
        end
        
        function writeHDF5Attributes(obj, filename, dataset)
            % Write attributes to HDF5 dataset
            
            % Write standard attributes
            h5writeatt(filename, dataset, 'created', datestr(now));
            h5writeatt(filename, dataset, 'creator', 'MUFO DataWriter');
            
            % Write metadata as attributes
            fields = fieldnames(obj.metadata);
            for i = 1:length(fields)
                value = obj.metadata.(fields{i});
                
                % Convert to writable format
                if isnumeric(value) || islogical(value)
                    h5writeatt(filename, dataset, fields{i}, value);
                elseif ischar(value) || isstring(value)
                    h5writeatt(filename, dataset, fields{i}, char(value));
                elseif isstruct(value)
                    % Flatten struct to attributes
                    obj.writeStructAsAttributes(filename, dataset, value, fields{i});
                end
            end
        end
        
        function writeStructAsAttributes(obj, filename, dataset, s, prefix)
            % Write structure fields as HDF5 attributes
            
            fields = fieldnames(s);
            for i = 1:length(fields)
                value = s.(fields{i});
                attrName = [prefix, '_', fields{i}];
                
                if isnumeric(value) || islogical(value) || ischar(value)
                    h5writeatt(filename, dataset, attrName, value);
                end
            end
        end
        
        function writeExchangeMetadata(obj, filename, metadata)
            % Write metadata in exchange format structure
            
            % Standard metadata locations
            if isfield(metadata, 'energy')
                obj.writeHDF5Dataset(filename, ...
                    '/measurement/instrument/source/energy', metadata.energy);
            end
            
            if isfield(metadata, 'pixelSize')
                obj.writeHDF5Dataset(filename, ...
                    '/measurement/instrument/detector/pixel_size', metadata.pixelSize);
            end
            
            if isfield(metadata, 'distance')
                obj.writeHDF5Dataset(filename, ...
                    '/measurement/instrument/detector/distance', metadata.distance);
            end
            
            if isfield(metadata, 'sampleName')
                obj.writeHDF5Dataset(filename, ...
                    '/measurement/sample/name', metadata.sampleName);
            end
        end
        
        function writeHDF5Dataset(obj, filename, dataset, data)
            % Write single value dataset
            
            % Create parent groups if needed
            parts = strsplit(dataset, '/');
            parts(cellfun(@isempty, parts)) = [];
            
            group = '/';
            for i = 1:length(parts)-1
                group = fullfile(group, parts{i});
                try
                    h5info(filename, group);
                catch
                    % Group doesn't exist, create it
                    if ~strcmp(group, '/')
                        % Get parent group
                        parentParts = strsplit(group, '/');
                        parentParts(cellfun(@isempty, parentParts)) = [];
                        if length(parentParts) > 1
                            parent = ['/', fullfile(parentParts{1:end-1})];
                        else
                            parent = '/';
                        end
                        
                        % Create group
                        plist = 'H5P_DEFAULT';
                        fid = H5F.open(filename, 'H5F_ACC_RDWR', plist);
                        gid = H5G.create(fid, group, plist, plist, plist);
                        H5G.close(gid);
                        H5F.close(fid);
                    end
                end
            end
            
            % Write dataset
            if ischar(data) || isstring(data)
                h5create(filename, dataset, 1, 'Datatype', 'string');
                h5write(filename, dataset, string(data));
            else
                h5create(filename, dataset, size(data));
                h5write(filename, dataset, data);
            end
        end
        
        function writeRAWMetadata(obj, metaFile, data, meta)
            % Write metadata for RAW file
            
            if nargin < 4
                % Extract from data
                meta = struct();
                meta.dims = size(data);
                meta.dataType = class(data);
                meta.byteOrder = obj.getSystemByteOrder();
            end
            
            % Write as JSON
            jsonStr = jsonencode(meta, 'PrettyPrint', true);
            
            fid = fopen(metaFile, 'w');
            fprintf(fid, '%s', jsonStr);
            fclose(fid);
        end
        
        function byteOrder = getSystemByteOrder(obj)
            % Get system byte order
            [~, ~, endian] = computer;
            if endian == 'L'
                byteOrder = 'little';
            else
                byteOrder = 'big';
            end
        end
        
        function merged = mergeStructs(obj, s1, s2)
            % Merge two structures
            merged = s1;
            fields = fieldnames(s2);
            
            for i = 1:length(fields)
                merged.(fields{i}) = s2.(fields{i});
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate DataWriter capabilities
            
            fprintf('\n=== DataWriter Demo ===\n\n');
            
            % Create writer
            writer = mufo.data.DataWriter('Compression', true);
            
            % Create test data
            testData = phantom(256) * 1000;
            testData = repmat(testData, [1, 1, 10]);
            
            % Example 1: Write HDF5
            fprintf('1. Writing HDF5 file:\n');
            writer.write(testData, 'test_output.h5', '/data');
            fprintf('   Written to test_output.h5\n');
            
            % Example 2: Write exchange format
            fprintf('\n2. Writing exchange format:\n');
            angles = linspace(0, 2*pi, 10);
            writer.writeExchangeFormat('test_exchange.h5', testData, ...
                'Angles', angles, ...
                'Flats', ones(256, 256, 5), ...
                'Darks', zeros(256, 256, 5));
            fprintf('   Written to test_exchange.h5\n');
            
            % Example 3: Incremental writing
            fprintf('\n3. Incremental writing:\n');
            writer.create('test_incremental.h5', [256, 256, 100], 'single');
            
            for i = 1:10
                chunk = rand(256, 256, 10, 'single');
                writer.writeSlices(chunk, (i-1)*10 + 1);
                fprintf('   Written chunk %d/10\n', i);
            end
            
            writer.close();
            
            % Clean up
            delete('test_output.h5');
            delete('test_exchange.h5');
            delete('test_incremental.h5');
            
            fprintf('\nDemo complete.\n');
        end
    end
end