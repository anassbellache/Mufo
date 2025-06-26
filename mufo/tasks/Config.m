classdef Config < handle
    % Config - JSON-based configuration management for MUFO
    %
    % This class handles loading, saving, and managing configuration
    % from JSON files, avoiding hardcoded paths and parameters.
    %
    % Example:
    %   cfg = mufo.Config('config.json');
    %   cfg.get('data.projections.path')
    %   cfg.set('reconstruction.center', 1024.5);
    %   cfg.save();
    
    properties (Access = private)
        configFile      % Path to JSON configuration file
        data            % Configuration data structure
        defaults        % Default configuration
        listeners       % Property change listeners
    end
    
    properties (SetAccess = private)
        isLoaded        % Whether config is loaded
        isModified      % Whether config has been modified
    end
    
    events
        ConfigLoaded
        ConfigModified
        ConfigSaved
    end
    
    methods
        function obj = Config(configFile)
            % Constructor
            %
            % Input:
            %   configFile - Path to JSON configuration file (optional)
            
            obj.isLoaded = false;
            obj.isModified = false;
            obj.listeners = {};
            
            % Set default configuration
            obj.defaults = obj.getDefaultConfig();
            obj.data = obj.defaults;
            
            % Load config file if provided
            if nargin > 0 && ~isempty(configFile)
                obj.load(configFile);
            end
        end
        
        function load(obj, configFile)
            % Load configuration from JSON file
            %
            % Input:
            %   configFile - Path to JSON file
            
            if ~exist(configFile, 'file')
                warning('Config:FileNotFound', ...
                    'Configuration file not found: %s\nCreating default configuration.', ...
                    configFile);
                obj.configFile = configFile;
                obj.save();
                return;
            end
            
            try
                % Read JSON file
                jsonStr = fileread(configFile);
                loadedData = jsondecode(jsonStr);
                
                % Merge with defaults
                obj.data = obj.mergeConfigs(obj.defaults, loadedData);
                
                obj.configFile = configFile;
                obj.isLoaded = true;
                obj.isModified = false;
                
                notify(obj, 'ConfigLoaded');
                
            catch ME
                error('Config:LoadError', ...
                    'Failed to load configuration: %s', ME.message);
            end
        end
        
        function save(obj, configFile)
            % Save configuration to JSON file
            %
            % Input:
            %   configFile - Path to save file (optional, uses loaded file)
            
            if nargin < 2
                configFile = obj.configFile;
            end
            
            if isempty(configFile)
                error('Config:NoFile', ...
                    'No configuration file specified');
            end
            
            try
                % Create directory if needed
                configDir = fileparts(configFile);
                if ~isempty(configDir) && ~exist(configDir, 'dir')
                    mkdir(configDir);
                end
                
                % Write JSON with pretty formatting
                jsonStr = jsonencode(obj.data, 'PrettyPrint', true);
                
                fid = fopen(configFile, 'w');
                fprintf(fid, '%s', jsonStr);
                fclose(fid);
                
                obj.configFile = configFile;
                obj.isModified = false;
                
                notify(obj, 'ConfigSaved');
                
            catch ME
                error('Config:SaveError', ...
                    'Failed to save configuration: %s', ME.message);
            end
        end
        
        function value = get(obj, path, default)
            % Get configuration value using dot notation
            %
            % Inputs:
            %   path - Dot-separated path (e.g., 'data.projections.path')
            %   default - Default value if path doesn't exist (optional)
            % Output:
            %   value - Configuration value
            
            if nargin < 3
                default = [];
            end
            
            parts = strsplit(path, '.');
            value = obj.data;
            
            try
                for i = 1:length(parts)
                    if isfield(value, parts{i})
                        value = value.(parts{i});
                    else
                        value = default;
                        return;
                    end
                end
            catch
                value = default;
            end
        end
        
        function set(obj, path, value)
            % Set configuration value using dot notation
            %
            % Inputs:
            %   path - Dot-separated path
            %   value - Value to set
            
            parts = strsplit(path, '.');
            
            % Navigate to parent
            current = obj.data;
            for i = 1:length(parts)-1
                if ~isfield(current, parts{i})
                    current.(parts{i}) = struct();
                end
                
                % Can't navigate further into non-struct
                if ~isstruct(current.(parts{i}))
                    error('Config:InvalidPath', ...
                        'Cannot set nested value in non-struct field: %s', ...
                        strjoin(parts(1:i), '.'));
                end
                
                if i < length(parts)-1
                    current = current.(parts{i});
                end
            end
            
            % Set value
            if length(parts) == 1
                obj.data.(parts{1}) = value;
            else
                parent = obj.getParent(parts);
                parent.(parts{end}) = value;
            end
            
            obj.isModified = true;
            notify(obj, 'ConfigModified');
        end
        
        function exists = has(obj, path)
            % Check if configuration path exists
            %
            % Input:
            %   path - Dot-separated path
            % Output:
            %   exists - Whether path exists
            
            parts = strsplit(path, '.');
            current = obj.data;
            
            for i = 1:length(parts)
                if isstruct(current) && isfield(current, parts{i})
                    current = current.(parts{i});
                else
                    exists = false;
                    return;
                end
            end
            
            exists = true;
        end
        
        function remove(obj, path)
            % Remove configuration value
            %
            % Input:
            %   path - Dot-separated path
            
            parts = strsplit(path, '.');
            
            if length(parts) == 1
                if isfield(obj.data, parts{1})
                    obj.data = rmfield(obj.data, parts{1});
                end
            else
                parent = obj.getParent(parts);
                if isfield(parent, parts{end})
                    parent = rmfield(parent, parts{end});
                    obj.setParent(parts, parent);
                end
            end
            
            obj.isModified = true;
            notify(obj, 'ConfigModified');
        end
        
        function config = getSection(obj, section)
            % Get entire configuration section
            %
            % Input:
            %   section - Section name (e.g., 'data', 'reconstruction')
            % Output:
            %   config - Section configuration
            
            config = obj.get(section, struct());
        end
        
        function setSection(obj, section, config)
            % Set entire configuration section
            %
            % Inputs:
            %   section - Section name
            %   config - Section configuration
            
            obj.set(section, config);
        end
        
        function paths = getPaths(obj)
            % Get all data paths from configuration
            paths = struct();
            
            % Extract common paths
            paths.projections = obj.get('data.projections.path', '');
            paths.darks = obj.get('data.darks.path', '');
            paths.flats = obj.get('data.flats.path', '');
            paths.output = obj.get('output.directory', '');
            paths.checkpoints = obj.get('output.checkpoints', '');
            paths.logs = obj.get('output.logs', '');
            
            % HDF5 dataset paths
            paths.h5Datasets = struct();
            paths.h5Datasets.projections = obj.get('data.projections.dataset', '/exchange/data');
            paths.h5Datasets.darks = obj.get('data.darks.dataset', '/exchange/data_dark');
            paths.h5Datasets.flats = obj.get('data.flats.dataset', '/exchange/data_white');
            paths.h5Datasets.theta = obj.get('data.angles.dataset', '/exchange/theta');
        end
        
        function params = getReconstructionParams(obj)
            % Get reconstruction parameters
            params = obj.getSection('reconstruction');
            
            % Ensure required fields exist
            if ~isfield(params, 'algorithm')
                params.algorithm = 'fbp';
            end
            if ~isfield(params, 'filter')
                params.filter = 'ramp-fromreal';
            end
            if ~isfield(params, 'center')
                params.center = [];  % Auto-detect
            end
        end
        
        function params = getPaganinParams(obj)
            % Get Paganin phase retrieval parameters
            params = obj.get('preprocessing.phaseRetrieval.paganin', struct());
            
            % Ensure required fields with defaults
            if ~isfield(params, 'energy')
                params.energy = 25;  % keV
            end
            if ~isfield(params, 'distance')
                params.distance = 1.0;  % meters
            end
            if ~isfield(params, 'pixelSize')
                params.pixelSize = 1e-6;  % meters
            end
            if ~isfield(params, 'delta')
                params.delta = 1e-7;
            end
            if ~isfield(params, 'beta')
                params.beta = 1e-9;
            end
            if ~isfield(params, 'applyLog')
                params.applyLog = true;
            end
        end
        
        function display(obj)
            % Display configuration in readable format
            fprintf('\n=== MUFO Configuration ===\n');
            fprintf('File: %s\n', obj.configFile);
            fprintf('Loaded: %s\n', mat2str(obj.isLoaded));
            fprintf('Modified: %s\n\n', mat2str(obj.isModified));
            
            obj.displayStruct(obj.data, 0);
        end
        
        function validate(obj)
            % Validate configuration
            errors = {};
            warnings = {};
            
            % Check required paths
            paths = obj.getPaths();
            
            if isempty(paths.projections)
                errors{end+1} = 'Missing projections path';
            elseif ~exist(paths.projections, 'file') && ~exist(paths.projections, 'dir')
                warnings{end+1} = sprintf('Projections path not found: %s', paths.projections);
            end
            
            if isempty(paths.output)
                errors{end+1} = 'Missing output directory';
            end
            
            % Check reconstruction parameters
            reconParams = obj.getReconstructionParams();
            
            validAlgorithms = {'fbp', 'sirt', 'cgls', 'sart', 'mlem', 'osem'};
            if ~any(strcmpi(reconParams.algorithm, validAlgorithms))
                warnings{end+1} = sprintf('Unknown reconstruction algorithm: %s', ...
                    reconParams.algorithm);
            end
            
            % Check Paganin parameters
            paganinParams = obj.getPaganinParams();
            
            if paganinParams.energy <= 0
                errors{end+1} = 'Paganin energy must be positive';
            end
            
            if paganinParams.distance <= 0
                errors{end+1} = 'Paganin distance must be positive';
            end
            
            if paganinParams.pixelSize <= 0
                errors{end+1} = 'Paganin pixel size must be positive';
            end
            
            % Report results
            if ~isempty(errors)
                fprintf('\n=== Configuration Errors ===\n');
                for i = 1:length(errors)
                    fprintf('ERROR: %s\n', errors{i});
                end
            end
            
            if ~isempty(warnings)
                fprintf('\n=== Configuration Warnings ===\n');
                for i = 1:length(warnings)
                    fprintf('WARNING: %s\n', warnings{i});
                end
            end
            
            if isempty(errors) && isempty(warnings)
                fprintf('Configuration is valid.\n');
            end
            
            % Return validation status
            if nargout > 0
                valid = isempty(errors);
            end
        end
        
        function template = createTemplate(obj, filename)
            % Create a template configuration file
            %
            % Input:
            %   filename - Output filename for template
            
            template = obj.getDefaultConfig();
            
            % Save template
            jsonStr = jsonencode(template, 'PrettyPrint', true);
            
            fid = fopen(filename, 'w');
            fprintf(fid, '%s', jsonStr);
            fclose(fid);
            
            fprintf('Configuration template created: %s\n', filename);
        end
    end
    
    methods (Access = private)
        function defaults = getDefaultConfig(obj)
            % Get default configuration structure
            defaults = struct();
            
            % Data paths
            defaults.data = struct();
            defaults.data.projections = struct('path', '', 'dataset', '/exchange/data');
            defaults.data.darks = struct('path', '', 'dataset', '/exchange/data_dark');
            defaults.data.flats = struct('path', '', 'dataset', '/exchange/data_white');
            defaults.data.angles = struct('dataset', '/exchange/theta', 'units', 'radians');
            
            % Preprocessing
            defaults.preprocessing = struct();
            defaults.preprocessing.flatField = struct(...
                'absorptionCorrect', true, ...
                'fixNanAndInf', true, ...
                'darkScale', 1.0, ...
                'flatScale', 1.0);
            
            defaults.preprocessing.phaseRetrieval = struct();
            defaults.preprocessing.phaseRetrieval.paganin = struct(...
                'enabled', false, ...
                'energy', 25, ...        % keV
                'distance', 1.0, ...     % meters
                'pixelSize', 1e-6, ...   % meters
                'delta', 1e-7, ...
                'beta', 1e-9, ...
                'applyLog', true);
            
            defaults.preprocessing.ringRemoval = struct(...
                'enabled', false, ...
                'method', 'wavelet', ...
                'parameters', struct());
            
            % Reconstruction
            defaults.reconstruction = struct();
            defaults.reconstruction.algorithm = 'fbp';
            defaults.reconstruction.filter = 'ramp-fromreal';
            defaults.reconstruction.center = [];  % Auto-detect
            defaults.reconstruction.angleRange = 360;
            defaults.reconstruction.precision = 'single';
            
            % Output
            defaults.output = struct();
            defaults.output.directory = 'output';
            defaults.output.format = 'tiff';
            defaults.output.dtype = 'float32';
            defaults.output.checkpoints = 'checkpoints';
            defaults.output.logs = 'logs';
            
            % UFO settings
            defaults.ufo = struct();
            defaults.ufo.executable = 'ufo-launch';
            defaults.ufo.pluginDir = '/usr/lib/x86_64-linux-gnu/ufo';
            defaults.ufo.kernelPath = '/usr/share/ufo/kernels';
            defaults.ufo.deviceType = 'gpu';
            defaults.ufo.numGpuThreads = [];  % Auto
            
            % Processing options
            defaults.processing = struct();
            defaults.processing.chunkSize = 500;
            defaults.processing.memoryLimit = [];  % Auto
            defaults.processing.verbose = true;
            defaults.processing.saveCheckpoints = false;
            defaults.processing.parallel = true;
            defaults.processing.numWorkers = [];  % Auto
        end
        
        function merged = mergeConfigs(obj, defaults, loaded)
            % Recursively merge loaded config with defaults
            merged = defaults;
            
            fields = fieldnames(loaded);
            for i = 1:length(fields)
                field = fields{i};
                
                if isstruct(defaults) && isfield(defaults, field) && ...
                   isstruct(defaults.(field)) && isstruct(loaded.(field))
                    % Recursively merge structs
                    merged.(field) = obj.mergeConfigs(defaults.(field), loaded.(field));
                else
                    % Overwrite with loaded value
                    merged.(field) = loaded.(field);
                end
            end
        end
        
        function parent = getParent(obj, parts)
            % Get parent struct for nested path
            parent = obj.data;
            
            for i = 1:length(parts)-1
                parent = parent.(parts{i});
            end
        end
        
        function setParent(obj, parts, value)
            % Set parent struct for nested path
            if length(parts) == 2
                obj.data.(parts{1}) = value;
            else
                parent = obj.data;
                for i = 1:length(parts)-2
                    parent = parent.(parts{i});
                end
                parent.(parts{end-1}) = value;
            end
        end
        
        function displayStruct(obj, s, indent)
            % Recursively display structure
            if nargin < 3
                indent = 0;
            end
            
            indentStr = repmat('  ', 1, indent);
            
            fields = fieldnames(s);
            for i = 1:length(fields)
                field = fields{i};
                value = s.(field);
                
                if isstruct(value)
                    fprintf('%s%s:\n', indentStr, field);
                    obj.displayStruct(value, indent + 1);
                elseif ischar(value) || isstring(value)
                    fprintf('%s%s: "%s"\n', indentStr, field, value);
                elseif isnumeric(value)
                    if isscalar(value)
                        fprintf('%s%s: %g\n', indentStr, field, value);
                    else
                        fprintf('%s%s: [%s]\n', indentStr, field, mat2str(value));
                    end
                elseif islogical(value)
                    fprintf('%s%s: %s\n', indentStr, field, mat2str(value));
                else
                    fprintf('%s%s: <%s>\n', indentStr, field, class(value));
                end
            end
        end
    end
end