classdef MetadataManager < handle
    % MetadataManager - Comprehensive metadata management for MUFO
    %
    % This class handles metadata throughout the tomographic reconstruction
    % pipeline, including acquisition parameters, processing history,
    % and data provenance tracking.
    %
    % Example:
    %   mm = mufo.data.MetadataManager();
    %   mm.loadFromFile('data.h5');
    %   mm.addProcessingStep('flat-field-correct', params);
    %   mm.saveToFile('processed.h5');
    
    properties (Access = private)
        metadata            % Main metadata structure
        schema              % Metadata schema for validation
        history             % Processing history stack
        fileInfo            % Source file information
        validationErrors    % Validation error log
    end
    
    properties (SetAccess = private)
        isModified          % Track if metadata has been modified
        schemaVersion       % Schema version number
    end
    
    properties (Constant)
        SCHEMA_VERSION = '1.0.0';
        
        % Standard metadata paths for HDF5 exchange format
        EXCHANGE_PATHS = struct(...
            'energy', '/measurement/instrument/source/energy', ...
            'pixelSize', '/measurement/instrument/detector/pixel_size', ...
            'distance', '/measurement/instrument/detector/distance', ...
            'exposure', '/measurement/instrument/detector/exposure_time', ...
            'sampleName', '/measurement/sample/name', ...
            'experimenter', '/measurement/experimenter/name', ...
            'facility', '/measurement/instrument/source/beamline', ...
            'date', '/measurement/date', ...
            'title', '/measurement/title', ...
            'angles', '/exchange/theta', ...
            'numAngles', '/process/acquisition/rotation/num_angles', ...
            'angleRange', '/process/acquisition/rotation/angle_range');
        
        % MATLAB toolbox metadata categories
        TOOLBOX_METADATA = struct(...
            'image', {'BitDepth', 'ColorType', 'Compression', 'Format'}, ...
            'dicom', {'StudyDate', 'Modality', 'PatientName', 'SliceThickness'}, ...
            'fits', {'SIMPLE', 'BITPIX', 'NAXIS', 'DATE-OBS'});
    end
    
    events
        MetadataLoaded
        MetadataModified
        MetadataSaved
        ValidationCompleted
        ProcessingStepAdded
    end
    
    methods
        function obj = MetadataManager(varargin)
            % Constructor
            %
            % Parameters:
            %   'Schema' - Custom schema file path
            %   'ValidateOnLoad' - Validate metadata on load (default: true)
            
            p = inputParser;
            addParameter(p, 'Schema', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'ValidateOnLoad', true, @islogical);
            parse(p, varargin{:});
            
            % Initialize
            obj.metadata = obj.createEmptyMetadata();
            obj.history = {};
            obj.validationErrors = {};
            obj.isModified = false;
            obj.schemaVersion = obj.SCHEMA_VERSION;
            
            % Load schema
            if ~isempty(p.Results.Schema)
                obj.loadSchema(p.Results.Schema);
            else
                obj.schema = obj.getDefaultSchema();
            end
        end
        
        %% Core Methods
        
        function loadFromFile(obj, filename, varargin)
            % Load metadata from file
            %
            % Inputs:
            %   filename - Source file path
            %   'Format' - Force specific format (auto-detect by default)
            
            p = inputParser;
            addRequired(p, 'filename', @(x) exist(x, 'file'));
            addParameter(p, 'Format', 'auto', @ischar);
            parse(p, filename, varargin{:});
            
            obj.fileInfo = dir(filename);
            [~, ~, ext] = fileparts(filename);
            
            % Auto-detect format
            format = p.Results.Format;
            if strcmpi(format, 'auto')
                format = obj.detectFormat(filename);
            end
            
            % Load based on format
            switch lower(format)
                case 'hdf5'
                    obj.loadHDF5Metadata(filename);
                case 'tiff'
                    obj.loadTIFFMetadata(filename);
                case 'dicom'
                    obj.loadDICOMMetadata(filename);
                case 'fits'
                    obj.loadFITSMetadata(filename);
                case 'json'
                    obj.loadJSONMetadata(filename);
                otherwise
                    warning('MetadataManager:UnknownFormat', ...
                        'Unknown format: %s', format);
            end
            
            obj.isModified = false;
            notify(obj, 'MetadataLoaded');
        end
        
        function saveToFile(obj, filename, varargin)
            % Save metadata to file
            %
            % Inputs:
            %   filename - Output file path
            %   'Format' - Output format (default: 'hdf5')
            %   'Overwrite' - Overwrite existing (default: false)
            
            p = inputParser;
            addRequired(p, 'filename');
            addParameter(p, 'Format', 'hdf5', @ischar);
            addParameter(p, 'Overwrite', false, @islogical);
            parse(p, filename, varargin{:});
            
            % Check if file exists
            if exist(filename, 'file') && ~p.Results.Overwrite
                error('MetadataManager:FileExists', ...
                    'File exists. Use ''Overwrite'', true to replace.');
            end
            
            % Save based on format
            switch lower(p.Results.Format)
                case 'hdf5'
                    obj.saveHDF5Metadata(filename);
                case 'json'
                    obj.saveJSONMetadata(filename);
                case 'xml'
                    obj.saveXMLMetadata(filename);
                otherwise
                    error('MetadataManager:UnsupportedFormat', ...
                        'Unsupported save format: %s', p.Results.Format);
            end
            
            obj.isModified = false;
            notify(obj, 'MetadataSaved');
        end
        
        function addProcessingStep(obj, taskName, parameters, varargin)
            % Add processing step to history
            %
            % Inputs:
            %   taskName - Name of the processing task
            %   parameters - Task parameters structure
            %   'Time' - Override timestamp
            %   'User' - User who performed the step
            %   'Notes' - Additional notes
            
            p = inputParser;
            addRequired(p, 'taskName', @ischar);
            addRequired(p, 'parameters', @isstruct);
            addParameter(p, 'Time', datetime('now'), @isdatetime);
            addParameter(p, 'User', getenv('USER'), @ischar);
            addParameter(p, 'Notes', '', @ischar);
            parse(p, taskName, parameters, varargin{:});
            
            % Create processing step entry
            step = struct();
            step.task = taskName;
            step.parameters = parameters;
            step.timestamp = p.Results.Time;
            step.user = p.Results.User;
            step.notes = p.Results.Notes;
            step.matlabVersion = version;
            step.platform = computer;
            
            % Add to history
            obj.history{end+1} = step;
            obj.metadata.processing.history = obj.history;
            obj.metadata.processing.lastModified = datetime('now');
            
            obj.isModified = true;
            notify(obj, 'ProcessingStepAdded');
            notify(obj, 'MetadataModified');
        end
        
        %% Getters and Setters
        
        function value = get(obj, key)
            % Get metadata value by key
            %
            % Input:
            %   key - Dot-separated key path (e.g., 'acquisition.energy')
            
            value = obj.getNestedField(obj.metadata, key);
        end
        
        function set(obj, key, value)
            % Set metadata value
            %
            % Inputs:
            %   key - Dot-separated key path
            %   value - Value to set
            
            obj.metadata = obj.setNestedField(obj.metadata, key, value);
            obj.isModified = true;
            notify(obj, 'MetadataModified');
        end
        
        function metadata = getMetadata(obj)
            % Get complete metadata structure
            metadata = obj.metadata;
        end
        
        function history = getProcessingHistory(obj)
            % Get processing history
            history = obj.history;
        end
        
        %% Validation Methods
        
        function [isValid, errors] = validate(obj, varargin)
            % Validate metadata against schema
            %
            % Parameters:
            %   'Strict' - Strict validation mode (default: false)
            %
            % Outputs:
            %   isValid - Validation result
            %   errors - List of validation errors
            
            p = inputParser;
            addParameter(p, 'Strict', false, @islogical);
            parse(p, varargin{:});
            
            errors = {};
            isValid = true;
            
            % Check required fields
            requiredFields = obj.schema.required;
            for i = 1:length(requiredFields)
                field = requiredFields{i};
                value = obj.get(field);
                
                if isempty(value)
                    errors{end+1} = sprintf('Required field missing: %s', field);
                    isValid = false;
                end
            end
            
            % Check field types
            fields = fieldnames(obj.schema.fields);
            for i = 1:length(fields)
                field = fields{i};
                value = obj.get(field);
                
                if ~isempty(value)
                    expectedType = obj.schema.fields.(field).type;
                    if ~obj.checkType(value, expectedType)
                        errors{end+1} = sprintf('Invalid type for %s: expected %s', ...
                            field, expectedType);
                        isValid = false;
                    end
                end
            end
            
            % Check ranges
            for i = 1:length(fields)
                field = fields{i};
                value = obj.get(field);
                fieldSchema = obj.schema.fields.(field);
                
                if ~isempty(value) && isfield(fieldSchema, 'range')
                    range = fieldSchema.range;
                    if isnumeric(value)
                        if value < range(1) || value > range(2)
                            errors{end+1} = sprintf('Value out of range for %s: %g not in [%g, %g]', ...
                                field, value, range(1), range(2));
                            isValid = false;
                        end
                    end
                end
            end
            
            obj.validationErrors = errors;
            notify(obj, 'ValidationCompleted');
        end
        
        %% Export Methods
        
        function exportReport(obj, filename, varargin)
            % Export metadata report
            %
            % Inputs:
            %   filename - Output file path
            %   'Format' - Report format: 'html', 'pdf', 'txt' (default: 'html')
            %   'IncludeHistory' - Include processing history (default: true)
            
            p = inputParser;
            addRequired(p, 'filename');
            addParameter(p, 'Format', 'html', @ischar);
            addParameter(p, 'IncludeHistory', true, @islogical);
            parse(p, filename, varargin{:});
            
            switch lower(p.Results.Format)
                case 'html'
                    obj.exportHTMLReport(filename, p.Results);
                case 'pdf'
                    obj.exportPDFReport(filename, p.Results);
                case 'txt'
                    obj.exportTextReport(filename, p.Results);
                otherwise
                    error('MetadataManager:UnsupportedFormat', ...
                        'Unsupported report format: %s', p.Results.Format);
            end
        end
        
        function summary = getSummary(obj)
            % Get metadata summary
            summary = struct();
            
            % Basic info
            summary.schemaVersion = obj.schemaVersion;
            summary.isModified = obj.isModified;
            summary.hasValidationErrors = ~isempty(obj.validationErrors);
            
            % Acquisition summary
            if isfield(obj.metadata, 'acquisition')
                acq = obj.metadata.acquisition;
                summary.acquisition = struct();
                
                if isfield(acq, 'energy')
                    summary.acquisition.energy = sprintf('%.1f keV', acq.energy);
                end
                
                if isfield(acq, 'numProjections')
                    summary.acquisition.projections = acq.numProjections;
                end
                
                if isfield(acq, 'pixelSize')
                    summary.acquisition.pixelSize = sprintf('%.2f µm', acq.pixelSize * 1e6);
                end
            end
            
            % Processing summary
            summary.processingSteps = length(obj.history);
            if ~isempty(obj.history)
                summary.lastProcessing = obj.history{end}.task;
                summary.lastProcessingTime = obj.history{end}.timestamp;
            end
        end
        
        %% Integration Methods
        
        function attachToDataLoader(obj, dataLoader)
            % Attach to DataLoader instance for automatic metadata extraction
            %
            % Input:
            %   dataLoader - mufo.data.DataLoader instance
            
            if ~isa(dataLoader, 'mufo.data.DataLoader')
                error('MetadataManager:InvalidInput', ...
                    'Input must be a DataLoader instance');
            end
            
            % Listen to DataLoader events
            addlistener(dataLoader, 'DataLoaded', ...
                @(src, evt) obj.extractFromDataLoader(src));
            addlistener(dataLoader, 'MetadataUpdated', ...
                @(src, evt) obj.mergeMetadata(src.getMetadata()));
        end
        
        function attachToDataWriter(obj, dataWriter)
            % Attach to DataWriter instance for automatic metadata writing
            %
            % Input:
            %   dataWriter - mufo.data.DataWriter instance
            
            if ~isa(dataWriter, 'mufo.data.DataWriter')
                error('MetadataManager:InvalidInput', ...
                    'Input must be a DataWriter instance');
            end
            
            % Set metadata in writer
            dataWriter.setMetadata(obj.getMetadata());
            
            % Listen to writer events
            addlistener(obj, 'MetadataModified', ...
                @(src, evt) dataWriter.setMetadata(obj.getMetadata()));
        end
        
        %% Private Methods
    end
    
    methods (Access = private)
        function metadata = createEmptyMetadata(obj)
            % Create empty metadata structure
            metadata = struct();
            
            % Standard categories
            metadata.acquisition = struct();
            metadata.instrument = struct();
            metadata.sample = struct();
            metadata.processing = struct();
            metadata.reconstruction = struct();
            
            % System info
            metadata.system = struct();
            metadata.system.created = datetime('now');
            metadata.system.matlabVersion = version;
            metadata.system.platform = computer;
            metadata.system.schemaVersion = obj.SCHEMA_VERSION;
        end
        
        function schema = getDefaultSchema(obj)
            % Get default metadata schema
            schema = struct();
            schema.version = obj.SCHEMA_VERSION;
            
            % Required fields
            schema.required = {'acquisition.numProjections', 'acquisition.pixelSize'};
            
            % Field definitions
            schema.fields = struct();
            
            % Acquisition fields
            schema.fields.energy = struct('type', 'numeric', 'range', [0.1, 1000], 'unit', 'keV');
            schema.fields.pixelSize = struct('type', 'numeric', 'range', [0.1e-6, 100e-6], 'unit', 'm');
            schema.fields.distance = struct('type', 'numeric', 'range', [0, 10], 'unit', 'm');
            schema.fields.exposure = struct('type', 'numeric', 'range', [0, 3600], 'unit', 's');
            schema.fields.numProjections = struct('type', 'integer', 'range', [1, 10000]);
            
            % Sample fields
            schema.fields.sampleName = struct('type', 'string', 'maxLength', 256);
            schema.fields.sampleDescription = struct('type', 'string', 'maxLength', 1024);
            
            % Processing fields
            schema.fields.centerOfRotation = struct('type', 'numeric', 'range', [0, inf]);
            schema.fields.ringRemovalMethod = struct('type', 'string', ...
                'enum', {'none', 'wavelet', 'morphological', 'combined'});
            
            return;
        end
        
        function loadHDF5Metadata(obj, filename)
            % Load metadata from HDF5 file
            
            info = h5info(filename);
            obj.metadata = obj.createEmptyMetadata();
            
            % Load exchange format metadata
            paths = fieldnames(obj.EXCHANGE_PATHS);
            for i = 1:length(paths)
                field = paths{i};
                h5path = obj.EXCHANGE_PATHS.(field);
                
                try
                    value = h5read(filename, h5path);
                    obj.set(['acquisition.', field], value);
                catch
                    % Path doesn't exist
                end
            end
            
            % Load attributes
            obj.loadHDF5Attributes(info, '');
            
            % Load custom metadata groups
            if any(strcmp({info.Groups.Name}, '/metadata'))
                metaGroup = h5info(filename, '/metadata');
                obj.loadHDF5Group(metaGroup, 'custom');
            end
        end
        
        function loadHDF5Attributes(obj, info, prefix)
            % Recursively load HDF5 attributes
            
            if isfield(info, 'Attributes')
                for i = 1:length(info.Attributes)
                    attr = info.Attributes(i);
                    key = [prefix, attr.Name];
                    obj.set(key, attr.Value);
                end
            end
            
            if isfield(info, 'Groups')
                for i = 1:length(info.Groups)
                    group = info.Groups(i);
                    groupName = strsplit(group.Name, '/');
                    groupName = groupName{end};
                    
                    newPrefix = [prefix, groupName, '.'];
                    obj.loadHDF5Attributes(group, newPrefix);
                end
            end
        end
        
        function loadHDF5Group(obj, groupInfo, prefix)
            % Load HDF5 group as metadata
            
            % Load datasets
            for i = 1:length(groupInfo.Datasets)
                dataset = groupInfo.Datasets(i);
                key = [prefix, '.', dataset.Name];
                
                try
                    value = h5read(groupInfo.Filename, ...
                        [groupInfo.Name, '/', dataset.Name]);
                    obj.set(key, value);
                catch
                    % Skip if can't read
                end
            end
            
            % Recurse into subgroups
            for i = 1:length(groupInfo.Groups)
                subgroup = groupInfo.Groups(i);
                subgroupName = strsplit(subgroup.Name, '/');
                subgroupName = subgroupName{end};
                
                obj.loadHDF5Group(subgroup, [prefix, '.', subgroupName]);
            end
        end
        
        function loadTIFFMetadata(obj, filename)
            % Load metadata from TIFF file using Image Processing Toolbox
            
            info = imfinfo(filename);
            obj.metadata = obj.createEmptyMetadata();
            
            % Basic TIFF metadata
            obj.set('format.type', 'TIFF');
            obj.set('format.width', info(1).Width);
            obj.set('format.height', info(1).Height);
            obj.set('format.bitDepth', info(1).BitDepth);
            obj.set('format.compression', info(1).Compression);
            
            if length(info) > 1
                obj.set('acquisition.numProjections', length(info));
            end
            
            % Extract TIFF tags
            if isfield(info(1), 'XResolution')
                pixelSize = 1 / info(1).XResolution;
                if isfield(info(1), 'ResolutionUnit') && ...
                   strcmpi(info(1).ResolutionUnit, 'Centimeter')
                    pixelSize = pixelSize * 0.01; % Convert to meters
                end
                obj.set('acquisition.pixelSize', pixelSize);
            end
            
            % Custom tags
            if isfield(info(1), 'UnknownTags')
                obj.parseTIFFTags(info(1).UnknownTags);
            end
            
            % ImageDescription field often contains metadata
            if isfield(info(1), 'ImageDescription')
                obj.parseImageDescription(info(1).ImageDescription);
            end
        end
        
        function loadDICOMMetadata(obj, filename)
            % Load metadata from DICOM file using Image Processing Toolbox
            
            if ~license('test', 'image_toolbox')
                warning('MetadataManager:NoToolbox', ...
                    'Image Processing Toolbox required for DICOM support');
                return;
            end
            
            info = dicominfo(filename);
            obj.metadata = obj.createEmptyMetadata();
            
            % Map DICOM fields to metadata
            dicomMap = struct(...
                'StudyDate', 'acquisition.date', ...
                'PatientName', 'sample.name', ...
                'KVP', 'acquisition.energy', ...
                'ExposureTime', 'acquisition.exposure', ...
                'PixelSpacing', 'acquisition.pixelSize', ...
                'SliceThickness', 'acquisition.sliceThickness');
            
            fields = fieldnames(dicomMap);
            for i = 1:length(fields)
                if isfield(info, fields{i})
                    value = info.(fields{i});
                    
                    % Convert specific fields
                    if strcmp(fields{i}, 'PixelSpacing') && isnumeric(value)
                        value = value(1) * 1e-3; % mm to m
                    elseif strcmp(fields{i}, 'ExposureTime')
                        value = value * 1e-3; % ms to s
                    end
                    
                    obj.set(dicomMap.(fields{i}), value);
                end
            end
        end
        
        function loadFITSMetadata(obj, filename)
            % Load metadata from FITS file
            
            info = fitsinfo(filename);
            obj.metadata = obj.createEmptyMetadata();
            
            % Process primary header
            if ~isempty(info.PrimaryData.Keywords)
                obj.parseFITSKeywords(info.PrimaryData.Keywords);
            end
            
            % Process image headers
            for i = 1:length(info.Image)
                if ~isempty(info.Image(i).Keywords)
                    prefix = sprintf('image%d.', i);
                    obj.parseFITSKeywords(info.Image(i).Keywords, prefix);
                end
            end
        end
        
        function parseFITSKeywords(obj, keywords, prefix)
            % Parse FITS keywords into metadata
            
            if nargin < 3
                prefix = '';
            end
            
            for i = 1:size(keywords, 1)
                key = keywords{i, 1};
                value = keywords{i, 2};
                
                % Skip standard FITS keywords
                if any(strcmpi(key, {'SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'COMMENT', 'HISTORY'}))
                    continue;
                end
                
                % Clean key name
                key = lower(strrep(key, '-', '_'));
                
                % Set value
                obj.set([prefix, key], value);
            end
        end
        
        function loadJSONMetadata(obj, filename)
            % Load metadata from JSON file
            
            fid = fopen(filename, 'r');
            jsonStr = fread(fid, '*char')';
            fclose(fid);
            
            obj.metadata = jsondecode(jsonStr);
            obj.isModified = false;
        end
        
        function saveHDF5Metadata(obj, filename)
            % Save metadata to HDF5 file
            
            % Check if file exists
            if exist(filename, 'file')
                % Append to existing file
                fileMode = 'append';
            else
                % Create new file
                fileMode = 'create';
            end
            
            % Write standard exchange format paths
            paths = fieldnames(obj.EXCHANGE_PATHS);
            for i = 1:length(paths)
                field = paths{i};
                value = obj.get(['acquisition.', field]);
                
                if ~isempty(value)
                    h5path = obj.EXCHANGE_PATHS.(field);
                    obj.writeHDF5Value(filename, h5path, value, fileMode);
                end
            end
            
            % Write custom metadata group
            obj.writeHDF5Group(filename, '/metadata', obj.metadata, fileMode);
            
            % Write processing history
            if ~isempty(obj.history)
                obj.writeProcessingHistory(filename);
            end
        end
        
        function writeHDF5Value(obj, filename, path, value, mode)
            % Write single value to HDF5
            
            try
                if strcmpi(mode, 'create') || ~obj.h5exists(filename, path)
                    % Create dataset
                    if ischar(value) || isstring(value)
                        h5create(filename, path, 1, 'Datatype', 'string');
                        h5write(filename, path, string(value));
                    else
                        h5create(filename, path, size(value));
                        h5write(filename, path, value);
                    end
                else
                    % Overwrite existing
                    h5write(filename, path, value);
                end
            catch ME
                warning('MetadataManager:WriteError', ...
                    'Failed to write %s: %s', path, ME.message);
            end
        end
        
        function writeHDF5Group(obj, filename, groupPath, data, mode)
            % Write structure as HDF5 group
            
            % Create group if needed
            if strcmpi(mode, 'create') || ~obj.h5exists(filename, groupPath)
                plist = 'H5P_DEFAULT';
                fid = H5F.open(filename, 'H5F_ACC_RDWR', plist);
                gid = H5G.create(fid, groupPath, plist, plist, plist);
                H5G.close(gid);
                H5F.close(fid);
            end
            
            % Write fields
            fields = fieldnames(data);
            for i = 1:length(fields)
                field = fields{i};
                value = data.(field);
                
                if isstruct(value)
                    % Recursive group creation
                    subgroupPath = [groupPath, '/', field];
                    obj.writeHDF5Group(filename, subgroupPath, value, mode);
                else
                    % Write dataset
                    datasetPath = [groupPath, '/', field];
                    obj.writeHDF5Value(filename, datasetPath, value, mode);
                end
            end
        end
        
        function writeProcessingHistory(obj, filename)
            % Write processing history to HDF5
            
            historyGroup = '/metadata/processing/history';
            
            for i = 1:length(obj.history)
                step = obj.history{i};
                stepGroup = sprintf('%s/step_%04d', historyGroup, i);
                
                % Write step data
                obj.writeHDF5Value(filename, [stepGroup, '/task'], step.task, 'append');
                obj.writeHDF5Value(filename, [stepGroup, '/timestamp'], ...
                    char(step.timestamp), 'append');
                obj.writeHDF5Value(filename, [stepGroup, '/user'], step.user, 'append');
                
                % Write parameters
                if ~isempty(fieldnames(step.parameters))
                    obj.writeHDF5Group(filename, [stepGroup, '/parameters'], ...
                        step.parameters, 'append');
                end
            end
        end
        
        function saveJSONMetadata(obj, filename)
            % Save metadata as JSON
            
            % Prepare data
            output = struct();
            output.metadata = obj.metadata;
            output.history = obj.history;
            output.schemaVersion = obj.schemaVersion;
            
            % Convert to JSON with pretty printing
            jsonStr = jsonencode(output, 'PrettyPrint', true);
            
            % Write file
            fid = fopen(filename, 'w');
            fprintf(fid, '%s', jsonStr);
            fclose(fid);
        end
        
        function saveXMLMetadata(obj, filename)
            % Save metadata as XML
            
            % Create DOM
            docNode = com.mathworks.xml.XMLUtils.createDocument('metadata');
            docRootNode = docNode.getDocumentElement;
            docRootNode.setAttribute('version', obj.schemaVersion);
            
            % Add metadata
            obj.struct2xml(docNode, docRootNode, obj.metadata);
            
            % Add history
            if ~isempty(obj.history)
                historyNode = docNode.createElement('history');
                docRootNode.appendChild(historyNode);
                
                for i = 1:length(obj.history)
                    stepNode = docNode.createElement('step');
                    stepNode.setAttribute('index', num2str(i));
                    obj.struct2xml(docNode, stepNode, obj.history{i});
                    historyNode.appendChild(stepNode);
                end
            end
            
            % Write file
            xmlwrite(filename, docNode);
        end
        
        function struct2xml(obj, docNode, parentNode, data)
            % Convert structure to XML nodes
            
            fields = fieldnames(data);
            for i = 1:length(fields)
                field = fields{i};
                value = data.(field);
                
                if isstruct(value)
                    % Create element and recurse
                    node = docNode.createElement(field);
                    parentNode.appendChild(node);
                    obj.struct2xml(docNode, node, value);
                else
                    % Create element with text
                    node = docNode.createElement(field);
                    
                    if isnumeric(value)
                        text = docNode.createTextNode(num2str(value));
                    elseif islogical(value)
                        text = docNode.createTextNode(mat2str(value));
                    elseif isdatetime(value)
                        text = docNode.createTextNode(char(value));
                    else
                        text = docNode.createTextNode(char(value));
                    end
                    
                    node.appendChild(text);
                    parentNode.appendChild(node);
                end
            end
        end
        
        function exportHTMLReport(obj, filename, options)
            % Export HTML report
            
            fid = fopen(filename, 'w');
            
            % HTML header
            fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n');
            fprintf(fid, '<title>Metadata Report</title>\n');
            fprintf(fid, '<style>\n');
            fprintf(fid, 'body { font-family: Arial, sans-serif; margin: 20px; }\n');
            fprintf(fid, 'table { border-collapse: collapse; width: 100%%; }\n');
            fprintf(fid, 'th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n');
            fprintf(fid, 'th { background-color: #4CAF50; color: white; }\n');
            fprintf(fid, 'h1, h2 { color: #333; }\n');
            fprintf(fid, '.error { color: red; }\n');
            fprintf(fid, '.warning { color: orange; }\n');
            fprintf(fid, '</style>\n</head>\n<body>\n');
            
            % Title
            fprintf(fid, '<h1>Tomography Metadata Report</h1>\n');
            fprintf(fid, '<p>Generated: %s</p>\n', char(datetime('now')));
            
            % Summary
            summary = obj.getSummary();
            fprintf(fid, '<h2>Summary</h2>\n');
            fprintf(fid, '<ul>\n');
            fprintf(fid, '<li>Schema Version: %s</li>\n', summary.schemaVersion);
            fprintf(fid, '<li>Processing Steps: %d</li>\n', summary.processingSteps);
            
            if isfield(summary, 'acquisition')
                fprintf(fid, '<li>Acquisition:</li><ul>\n');
                fields = fieldnames(summary.acquisition);
                for i = 1:length(fields)
                    fprintf(fid, '<li>%s: %s</li>\n', fields{i}, ...
                        string(summary.acquisition.(fields{i})));
                end
                fprintf(fid, '</ul>\n');
            end
            fprintf(fid, '</ul>\n');
            
            % Validation
            if ~isempty(obj.validationErrors)
                fprintf(fid, '<h2>Validation Errors</h2>\n');
                fprintf(fid, '<ul class="error">\n');
                for i = 1:length(obj.validationErrors)
                    fprintf(fid, '<li>%s</li>\n', obj.validationErrors{i});
                end
                fprintf(fid, '</ul>\n');
            end
            
            % Metadata table
            fprintf(fid, '<h2>Metadata</h2>\n');
            fprintf(fid, '<table>\n<tr><th>Key</th><th>Value</th></tr>\n');
            obj.writeMetadataTable(fid, obj.metadata, '');
            fprintf(fid, '</table>\n');
            
            % Processing history
            if options.IncludeHistory && ~isempty(obj.history)
                fprintf(fid, '<h2>Processing History</h2>\n');
                fprintf(fid, '<table>\n');
                fprintf(fid, '<tr><th>Step</th><th>Task</th><th>Time</th><th>User</th></tr>\n');
                
                for i = 1:length(obj.history)
                    step = obj.history{i};
                    fprintf(fid, '<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td></tr>\n', ...
                        i, step.task, char(step.timestamp), step.user);
                end
                
                fprintf(fid, '</table>\n');
            end
            
            % Footer
            fprintf(fid, '</body>\n</html>\n');
            fclose(fid);
        end
        
        function writeMetadataTable(obj, fid, data, prefix)
            % Write metadata as HTML table rows
            
            fields = fieldnames(data);
            for i = 1:length(fields)
                field = fields{i};
                value = data.(field);
                key = [prefix, field];
                
                if isstruct(value)
                    % Recurse
                    obj.writeMetadataTable(fid, value, [key, '.']);
                else
                    % Write row
                    if isnumeric(value)
                        valueStr = num2str(value);
                    elseif islogical(value)
                        valueStr = mat2str(value);
                    else
                        valueStr = char(value);
                    end
                    
                    fprintf(fid, '<tr><td>%s</td><td>%s</td></tr>\n', key, valueStr);
                end
            end
        end
        
        function exportTextReport(obj, filename, options)
            % Export text report
            
            fid = fopen(filename, 'w');
            
            fprintf(fid, 'TOMOGRAPHY METADATA REPORT\n');
            fprintf(fid, '=========================\n\n');
            fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));
            
            % Summary
            summary = obj.getSummary();
            fprintf(fid, 'SUMMARY\n-------\n');
            fprintf(fid, 'Schema Version: %s\n', summary.schemaVersion);
            fprintf(fid, 'Processing Steps: %d\n', summary.processingSteps);
            
            % Metadata
            fprintf(fid, '\nMETADATA\n--------\n');
            obj.writeMetadataText(fid, obj.metadata, '');
            
            % History
            if options.IncludeHistory && ~isempty(obj.history)
                fprintf(fid, '\nPROCESSING HISTORY\n-----------------\n');
                for i = 1:length(obj.history)
                    step = obj.history{i};
                    fprintf(fid, '\nStep %d: %s\n', i, step.task);
                    fprintf(fid, '  Time: %s\n', char(step.timestamp));
                    fprintf(fid, '  User: %s\n', step.user);
                end
            end
            
            fclose(fid);
        end
        
        function writeMetadataText(obj, fid, data, prefix)
            % Write metadata as text
            
            fields = fieldnames(data);
            for i = 1:length(fields)
                field = fields{i};
                value = data.(field);
                key = [prefix, field];
                
                if isstruct(value)
                    fprintf(fid, '\n%s:\n', key);
                    obj.writeMetadataText(fid, value, [prefix, '  ']);
                else
                    if isnumeric(value)
                        fprintf(fid, '%s: %g\n', key, value);
                    else
                        fprintf(fid, '%s: %s\n', key, char(value));
                    end
                end
            end
        end
        
        function value = getNestedField(obj, data, key)
            % Get nested field value using dot notation
            
            parts = strsplit(key, '.');
            value = data;
            
            for i = 1:length(parts)
                if isstruct(value) && isfield(value, parts{i})
                    value = value.(parts{i});
                else
                    value = [];
                    return;
                end
            end
        end
        
        function data = setNestedField(obj, data, key, value)
            % Set nested field value using dot notation
            
            parts = strsplit(key, '.');
            
            % Build nested structure
            if length(parts) == 1
                data.(parts{1}) = value;
            else
                if ~isfield(data, parts{1})
                    data.(parts{1}) = struct();
                end
                
                subkey = strjoin(parts(2:end), '.');
                data.(parts{1}) = obj.setNestedField(data.(parts{1}), subkey, value);
            end
        end
        
        function isValid = checkType(obj, value, expectedType)
            % Check if value matches expected type
            
            switch expectedType
                case 'numeric'
                    isValid = isnumeric(value);
                case 'integer'
                    isValid = isnumeric(value) && all(value == round(value));
                case 'string'
                    isValid = ischar(value) || isstring(value);
                case 'logical'
                    isValid = islogical(value);
                case 'datetime'
                    isValid = isdatetime(value);
                otherwise
                    isValid = true;
            end
        end
        
        function format = detectFormat(obj, filename)
            % Auto-detect file format
            
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case {'.h5', '.hdf5', '.nxs', '.nx'}
                    format = 'hdf5';
                case {'.tif', '.tiff'}
                    format = 'tiff';
                case {'.dcm', '.dicom'}
                    format = 'dicom';
                case {'.fits', '.fit'}
                    format = 'fits';
                case '.json'
                    format = 'json';
                case '.xml'
                    format = 'xml';
                otherwise
                    % Try to detect by content
                    try
                        h5info(filename);
                        format = 'hdf5';
                    catch
                        format = 'unknown';
                    end
            end
        end
        
        function exists = h5exists(obj, filename, path)
            % Check if HDF5 path exists
            
            try
                h5info(filename, path);
                exists = true;
            catch
                exists = false;
            end
        end
        
        function extractFromDataLoader(obj, dataLoader)
            % Extract metadata from DataLoader
            
            loaderMeta = dataLoader.getMetadata();
            
            % Merge with existing metadata
            obj.mergeMetadata(loaderMeta);
        end
        
        function mergeMetadata(obj, newMetadata)
            % Merge new metadata with existing
            
            obj.metadata = obj.mergeStructs(obj.metadata, newMetadata);
            obj.isModified = true;
            notify(obj, 'MetadataModified');
        end
        
        function merged = mergeStructs(obj, s1, s2)
            % Recursively merge two structures
            
            merged = s1;
            fields = fieldnames(s2);
            
            for i = 1:length(fields)
                field = fields{i};
                
                if isfield(merged, field) && isstruct(merged.(field)) && isstruct(s2.(field))
                    % Recursive merge
                    merged.(field) = obj.mergeStructs(merged.(field), s2.(field));
                else
                    % Direct assignment
                    merged.(field) = s2.(field);
                end
            end
        end
    end
    
    %% Static Methods
    methods (Static)
        function demo()
            % Demonstrate MetadataManager capabilities
            
            fprintf('\n=== MetadataManager Demo ===\n\n');
            
            % Create manager
            mm = mufo.data.MetadataManager();
            
            % Set basic metadata
            fprintf('1. Setting metadata:\n');
            mm.set('acquisition.energy', 25.0);
            mm.set('acquisition.pixelSize', 1.1e-6);
            mm.set('acquisition.numProjections', 1800);
            mm.set('sample.name', 'Test Sample');
            
            summary = mm.getSummary();
            fprintf('   Energy: %s\n', summary.acquisition.energy);
            fprintf('   Pixel Size: %s\n', summary.acquisition.pixelSize);
            
            % Add processing step
            fprintf('\n2. Adding processing step:\n');
            params = struct('method', 'wavelet', 'level', 3);
            mm.addProcessingStep('ring-remove', params);
            fprintf('   Processing steps: %d\n', length(mm.getProcessingHistory()));
            
            % Validate
            fprintf('\n3. Validating metadata:\n');
            [isValid, errors] = mm.validate();
            if isValid
                fprintf('   Validation passed!\n');
            else
                fprintf('   Validation errors:\n');
                for i = 1:length(errors)
                    fprintf('     - %s\n', errors{i});
                end
            end
            
            % Save to different formats
            fprintf('\n4. Saving metadata:\n');
            mm.saveToFile('test_metadata.json', 'Format', 'json', 'Overwrite', true);
            fprintf('   Saved to test_metadata.json\n');
            
            % Export report
            fprintf('\n5. Exporting report:\n');
            mm.exportReport('test_report.html', 'Format', 'html');
            fprintf('   Report saved to test_report.html\n');
            
            % Load from file
            fprintf('\n6. Loading metadata:\n');
            mm2 = mufo.data.MetadataManager();
            mm2.loadFromFile('test_metadata.json');
            fprintf('   Loaded %d processing steps\n', length(mm2.getProcessingHistory()));
            
            % Clean up
            delete('test_metadata.json');
            delete('test_report.html');
            
            fprintf('\nDemo complete!\n');
        end
    end
end