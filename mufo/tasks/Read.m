classdef Read < mufo.core.UFOTask
    % Read - UFO read task wrapper with HDF5 support
    %
    % This task reads data from various formats including HDF5, TIFF,
    % and raw files. It integrates with MATLAB's HDF5 capabilities
    % for dataset inspection and selection.
    %
    % Example:
    %   task = mufo.tasks.Read();
    %   task.setPath('data.h5:/exchange/data');
    %   task.setRegion(100, 200, 50, 100);  % x, y, width, height
    
    properties (Access = private)
        dataInfo        % Information about the data source
    end
    
    methods
        function obj = Read(varargin)
            % Constructor
            obj@mufo.core.UFOTask('read');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setPath(obj, path)
            % Set data path with automatic format detection
            %
            % Input:
            %   path - File path or pattern (supports HDF5 notation)
            
            % Parse HDF5 notation (file.h5:/dataset/path)
            if contains(path, ':/')
                parts = strsplit(path, ':/');
                filepath = parts{1};
                dataset = parts{2};
                
                % Validate HDF5 file
                if exist(filepath, 'file')
                    obj.setParameter('path', filepath);
                    obj.setParameter('hdf5-dataset', dataset);
                    
                    % Get dataset info
                    obj.dataInfo = obj.inspectHDF5(filepath, dataset);
                else
                    error('Read:FileNotFound', 'HDF5 file not found: %s', filepath);
                end
            else
                % Regular file or pattern
                obj.setParameter('path', path);
                
                % Auto-detect format
                [~, ~, ext] = fileparts(path);
                if strcmpi(ext, '.h5') || strcmpi(ext, '.hdf5')
                    % Try to find main dataset
                    dataset = obj.findMainDataset(path);
                    if ~isempty(dataset)
                        obj.setParameter('hdf5-dataset', dataset);
                        obj.dataInfo = obj.inspectHDF5(path, dataset);
                    end
                end
            end
        end
        
        function setRegion(obj, x, y, width, height)
            % Set region of interest to read
            %
            % Inputs:
            %   x, y - Starting position (0-based for UFO)
            %   width, height - Region size
            
            obj.setParameter('x', x);
            obj.setParameter('y', y);
            obj.setParameter('width', width);
            obj.setParameter('height', height);
        end
        
        function setSliceRange(obj, start, count, step)
            % Set slice range for 3D data
            %
            % Inputs:
            %   start - Starting slice (0-based)
            %   count - Number of slices
            %   step - Step between slices (default: 1)
            
            obj.setParameter('start', start);
            obj.setParameter('number', count);
            
            if nargin > 3 && step > 1
                obj.setParameter('step', step);
            end
        end
        
        function info = inspectData(obj)
            % Inspect data source and return information
            %
            % Output:
            %   info - Structure with data information
            
            path = obj.getParameter('path');
            
            if isempty(path)
                error('Read:NoPath', 'No data path set');
            end
            
            % Check if it's HDF5
            [~, ~, ext] = fileparts(path);
            if strcmpi(ext, '.h5') || strcmpi(ext, '.hdf5')
                dataset = obj.getParameter('hdf5-dataset');
                if isempty(dataset)
                    dataset = obj.findMainDataset(path);
                end
                info = obj.inspectHDF5(path, dataset);
            else
                % Try to get info from file pattern
                info = obj.inspectFilePattern(path);
            end
            
            obj.dataInfo = info;
        end
        
        function preview = previewData(obj, sliceIdx)
            % Preview data using MATLAB's visualization
            %
            % Input:
            %   sliceIdx - Slice index to preview (optional)
            % Output:
            %   preview - Preview data
            
            if isempty(obj.dataInfo)
                obj.inspectData();
            end
            
            path = obj.getParameter('path');
            
            % Read preview data
            if obj.dataInfo.isHDF5
                dataset = obj.getParameter('hdf5-dataset');
                
                if nargin < 2
                    sliceIdx = round(obj.dataInfo.dims(3) / 2);
                end
                
                % Read single slice using MATLAB's HDF5 functions
                start = [1, 1, sliceIdx];
                count = [obj.dataInfo.dims(1), obj.dataInfo.dims(2), 1];
                preview = h5read(path, dataset, start, count);
                preview = squeeze(preview);
                
                % Display
                figure('Name', sprintf('Data Preview - Slice %d', sliceIdx));
                imagesc(preview);
                colormap gray;
                axis image;
                colorbar;
                title(sprintf('%s - Slice %d/%d', dataset, sliceIdx, obj.dataInfo.dims(3)));
                
            else
                % Read first image from pattern
                files = dir(path);
                if ~isempty(files)
                    preview = imread(fullfile(files(1).folder, files(1).name));
                    
                    figure('Name', 'Data Preview');
                    imagesc(preview);
                    colormap gray;
                    axis image;
                    colorbar;
                    title(sprintf('First image: %s', files(1).name));
                else
                    error('Read:NoFiles', 'No files found matching pattern');
                end
            end
        end
        
        function datasets = listHDF5Datasets(obj, filepath)
            % List all datasets in HDF5 file
            %
            % Input:
            %   filepath - HDF5 file path
            % Output:
            %   datasets - Cell array of dataset paths
            
            if nargin < 2
                filepath = obj.getParameter('path');
            end
            
            if ~exist(filepath, 'file')
                error('Read:FileNotFound', 'File not found: %s', filepath);
            end
            
            info = h5info(filepath);
            datasets = {};
            
            % Recursively find datasets
            datasets = obj.findDatasetsRecursive(info, '', datasets);
            
            % Display if no output requested
            if nargout == 0
                fprintf('\nDatasets in %s:\n', filepath);
                for i = 1:length(datasets)
                    dsInfo = h5info(filepath, datasets{i});
                    fprintf('  %s [%s]\n', datasets{i}, ...
                        strjoin(arrayfun(@num2str, dsInfo.Dataspace.Size, 'UniformOutput', false), 'x'));
                end
            end
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                paths = config.getPaths();
                
                % Set projection path
                if ~isempty(paths.projections)
                    if contains(paths.projections, '.h5')
                        fullPath = sprintf('%s:%s', paths.projections, ...
                            paths.h5Datasets.projections);
                        obj.setPath(fullPath);
                    else
                        obj.setPath(paths.projections);
                    end
                end
            else
                % Direct struct configuration
                if isfield(config, 'path')
                    obj.setPath(config.path);
                end
                
                if isfield(config, 'region')
                    obj.setRegion(config.region.x, config.region.y, ...
                        config.region.width, config.region.height);
                end
                
                if isfield(config, 'slices')
                    obj.setSliceRange(config.slices.start, ...
                        config.slices.count, config.slices.step);
                end
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define read task parameters
            
            % Path (required)
            obj.addParameterDef('path', 'string', ...
                'description', 'File path or pattern', ...
                'required', true);
            
            % HDF5 specific
            obj.addParameterDef('hdf5-dataset', 'string', ...
                'description', 'HDF5 dataset path', ...
                'default', '');
            
            % Region of interest
            obj.addParameterDef('x', 'integer', ...
                'description', 'X offset', ...
                'default', 0, ...
                'range', [0, inf]);
            
            obj.addParameterDef('y', 'integer', ...
                'description', 'Y offset', ...
                'default', 0, ...
                'range', [0, inf]);
            
            obj.addParameterDef('width', 'integer', ...
                'description', 'Width to read', ...
                'range', [1, inf]);
            
            obj.addParameterDef('height', 'integer', ...
                'description', 'Height to read', ...
                'range', [1, inf]);
            
            % Slice selection
            obj.addParameterDef('start', 'integer', ...
                'description', 'Starting slice', ...
                'default', 0, ...
                'range', [0, inf]);
            
            obj.addParameterDef('number', 'integer', ...
                'description', 'Number of items to read', ...
                'range', [1, inf]);
            
            obj.addParameterDef('step', 'integer', ...
                'description', 'Step between items', ...
                'default', 1, ...
                'range', [1, inf]);
            
            % Data type conversion
            obj.addParameterDef('convert', 'string', ...
                'description', 'Convert data type', ...
                'values', {'', 'float', 'uint8', 'uint16'});
            
            % Byte order
            obj.addParameterDef('swab', 'logical', ...
                'description', 'Swap byte order', ...
                'default', false);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            % Check if path is set
            if isempty(obj.getParameter('path'))
                valid = false;
                return;
            end
            
            % Validate region parameters
            x = obj.getParameter('x');
            y = obj.getParameter('y');
            width = obj.getParameter('width');
            height = obj.getParameter('height');
            
            if ~isempty(width) && ~isempty(height)
                if width <= 0 || height <= 0
                    valid = false;
                    warning('Read:InvalidRegion', 'Width and height must be positive');
                end
            end
        end
    end
    
    methods (Access = private)
        function info = inspectHDF5(obj, filepath, dataset)
            % Inspect HDF5 dataset
            info = struct();
            info.isHDF5 = true;
            info.file = filepath;
            info.dataset = dataset;
            
            try
                h5info_data = h5info(filepath, dataset);
                info.dims = h5info_data.Dataspace.Size;
                info.dtype = h5info_data.Datatype.Class;
                
                % Get attributes
                info.attributes = struct();
                for i = 1:length(h5info_data.Attributes)
                    attr = h5info_data.Attributes(i);
                    info.attributes.(attr.Name) = attr.Value;
                end
                
                % Check for common metadata
                try
                    info.theta = h5read(filepath, '/exchange/theta');
                catch
                    info.theta = [];
                end
                
            catch ME
                warning('Read:HDF5Error', 'Error inspecting HDF5: %s', ME.message);
                info = struct('isHDF5', true, 'error', ME.message);
            end
        end
        
        function info = inspectFilePattern(obj, pattern)
            % Inspect file pattern
            info = struct();
            info.isHDF5 = false;
            info.pattern = pattern;
            
            % Find matching files
            files = dir(pattern);
            info.numFiles = length(files);
            
            if info.numFiles > 0
                % Read first file to get dimensions
                firstFile = fullfile(files(1).folder, files(1).name);
                try
                    img = imread(firstFile);
                    info.dims = [size(img, 2), size(img, 1), info.numFiles];
                    info.dtype = class(img);
                catch
                    info.dims = [0, 0, info.numFiles];
                    info.dtype = 'unknown';
                end
            else
                info.dims = [0, 0, 0];
                info.dtype = 'unknown';
            end
        end
        
        function dataset = findMainDataset(obj, filepath)
            % Find main dataset in HDF5 file
            datasets = obj.listHDF5Datasets(filepath);
            
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
                    info = h5info(filepath, datasets{i});
                    dataSize = prod(info.Dataspace.Size);
                    if dataSize > maxSize
                        maxSize = dataSize;
                        dataset = datasets{i};
                    end
                end
            else
                dataset = '';
            end
        end
        
        function datasets = findDatasetsRecursive(obj, info, prefix, datasets)
            % Recursively find all datasets in HDF5 structure
            
            % Check datasets at current level
            for i = 1:length(info.Datasets)
                datasetPath = fullfile(prefix, info.Datasets(i).Name);
                datasets{end+1} = datasetPath;
            end
            
            % Recursively check groups
            for i = 1:length(info.Groups)
                groupPath = info.Groups(i).Name;
                datasets = obj.findDatasetsRecursive(info.Groups(i), ...
                    groupPath, datasets);
            end
        end
    end
end