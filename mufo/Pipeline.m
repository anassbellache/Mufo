classdef Pipeline < mufo.core.UFOChain
    % Pipeline - High-level pipeline orchestrator for MUFO
    %
    % This class extends UFOChain with additional features for
    % interactive inspection, checkpoint management, and integration
    % with MATLAB visualization tools.
    %
    % Example:
    %   pipeline = mufo.Pipeline(config);
    %   pipeline.addTask(task, 'checkpoint', true);
    %   pipeline.run('interactive', true);
    
    properties (Access = private)
        config          % Configuration object
        cache           % Data cache manager
        visualizer      % Visualization engine
        logger          % Logging system
        checkpoints     % Checkpoint information
        tempFiles       % Temporary file tracking
    end
    
    properties
        enableCache = true      % Enable data caching
        enableLogging = true    % Enable logging
        cleanupTemp = true      % Cleanup temporary files
    end
    
    events
        CheckpointReached
        VisualizationRequested
        PipelineCompleted
    end
    
    methods
        function obj = Pipeline(config, varargin)
            % Constructor
            %
            % Inputs:
            %   config - Config object or path to config file
            %   varargin - Additional options for UFOChain
            
            % Initialize parent
            obj@mufo.core.UFOChain(varargin{:});
            
            % Load configuration
            if ischar(config) || isstring(config)
                obj.config = mufo.Config(config);
            elseif isa(config, 'mufo.Config')
                obj.config = config;
            else
                error('Pipeline:InvalidConfig', ...
                    'Config must be a Config object or path to config file');
            end
            
            % Initialize components
            obj.initializeComponents();
            
            % Set UFO environment from config
            obj.setupEnvironment();
        end
        
        function addTask(obj, task, varargin)
            % Override to add checkpoint support
            %
            % Additional parameters:
            %   'checkpoint' - Add checkpoint after task
            %   'name' - Checkpoint name
            %   'visualize' - Visualization function
            
            p = inputParser;
            addRequired(p, 'task');
            addParameter(p, 'checkpoint', false, @islogical);
            addParameter(p, 'name', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'visualize', [], @(x) isa(x, 'function_handle'));
            
            % Parse known parameters
            params = varargin;
            checkpointIdx = find(strcmpi(params, 'checkpoint'));
            if ~isempty(checkpointIdx)
                checkpoint = params{checkpointIdx + 1};
                params([checkpointIdx, checkpointIdx + 1]) = [];
            else
                checkpoint = false;
            end
            
            nameIdx = find(strcmpi(params, 'name'));
            if ~isempty(nameIdx)
                name = params{nameIdx + 1};
                params([nameIdx, nameIdx + 1]) = [];
            else
                name = '';
            end
            
            visualizeIdx = find(strcmpi(params, 'visualize'));
            if ~isempty(visualizeIdx)
                visualizeFcn = params{visualizeIdx + 1};
                params([visualizeIdx, visualizeIdx + 1]) = [];
            else
                visualizeFcn = [];
            end
            
            % Add task to chain
            addTask@mufo.core.UFOChain(obj, task, params{:});
            
            % Store checkpoint info
            if checkpoint
                obj.checkpoints{end+1} = struct(...
                    'index', obj.length, ...
                    'name', name, ...
                    'visualize', visualizeFcn);
            end
        end
        
        function [status, output] = run(obj, varargin)
            % Execute pipeline with enhanced features
            %
            % Parameters:
            %   'interactive' - Enable interactive inspection
            %   'parallel' - Use parallel processing where possible
            %   'profile' - Enable performance profiling
            
            p = inputParser;
            addParameter(p, 'interactive', false, @islogical);
            addParameter(p, 'parallel', true, @islogical);
            addParameter(p, 'profile', false, @islogical);
            parse(p, varargin{:});
            
            % Start logging
            if obj.enableLogging
                obj.logger.startLog('pipeline_execution');
            end
            
            % Start profiling
            if p.Results.profile
                profile on;
            end
            
            try
                if p.Results.interactive && ~isempty(obj.checkpoints)
                    % Interactive execution with checkpoints
                    [status, output] = obj.runInteractive();
                else
                    % Standard execution
                    [status, output] = obj.execute('Process', true);
                end
                
                notify(obj, 'PipelineCompleted');
                
            catch ME
                obj.logger.error('Pipeline execution failed: %s', ME.message);
                rethrow(ME);
            end
            
            % Stop profiling
            if p.Results.profile
                profile off;
                profData = profile('info');
                obj.saveProfilingReport(profData);
            end
            
            % Cleanup
            if obj.cleanupTemp
                obj.cleanupTempFiles();
            end
        end
        
        function data = loadCheckpoint(obj, checkpointName)
            % Load data from checkpoint
            %
            % Input:
            %   checkpointName - Name of checkpoint
            % Output:
            %   data - Loaded data
            
            if obj.enableCache
                data = obj.cache.load(checkpointName);
            else
                error('Pipeline:CacheDisabled', 'Cache is disabled');
            end
        end
        
        function visualizeData(obj, data, options)
            % Visualize data using appropriate viewer
            %
            % Inputs:
            %   data - Data to visualize
            %   options - Visualization options
            
            if nargin < 3
                options = struct();
            end
            
            % Auto-detect data type
            if obj.isSinogram(data)
                obj.visualizer.showSinogram(data, options);
            elseif obj.isProjection(data)
                obj.visualizer.showProjection(data, options);
            elseif obj.isVolume(data)
                obj.visualizer.showVolume(data, options);
            else
                obj.visualizer.showGeneric(data, options);
            end
            
            notify(obj, 'VisualizationRequested');
        end
        
        function report = generateReport(obj)
            % Generate execution report
            %
            % Output:
            %   report - Report structure
            
            report = struct();
            report.pipeline = obj.getChainInfo();
            report.config = struct(obj.config.data);
            report.execution = obj.logger.getExecutionStats();
            
            if obj.enableCache
                report.cache = obj.cache.getStats();
            end
            
            % Save report
            reportFile = fullfile(obj.config.get('output.logs'), ...
                sprintf('pipeline_report_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
            save(reportFile, 'report');
            
            fprintf('Report saved to: %s\n', reportFile);
        end
    end
    
    methods (Access = private)
        function initializeComponents(obj)
            % Initialize pipeline components
            
            % Data cache
            if obj.enableCache
                cacheDir = obj.config.get('output.checkpoints', tempdir);
                obj.cache = mufo.data.DataCache(cacheDir);
            end
            
            % Visualizer
            obj.visualizer = mufo.vis.Inspector();
            
            % Logger
            if obj.enableLogging
                logDir = obj.config.get('output.logs', tempdir);
                obj.logger = mufo.Logger('LogDir', logDir);
            end
            
            % Initialize checkpoint list
            obj.checkpoints = {};
            obj.tempFiles = {};
        end
        
        function setupEnvironment(obj)
            % Setup UFO environment from config
            
            ufoConfig = obj.config.getSection('ufo');
            
            if isfield(ufoConfig, 'pluginDir')
                obj.command.setEnvironment('UFO_PLUGIN_DIR', ufoConfig.pluginDir);
            end
            
            if isfield(ufoConfig, 'kernelPath')
                obj.command.setEnvironment('UFO_KERNEL_PATH', ufoConfig.kernelPath);
            end
            
            if isfield(ufoConfig, 'deviceType')
                obj.command.setEnvironment('UFO_DEVICE_TYPE', ufoConfig.deviceType);
            end
        end
        
        function [status, output] = runInteractive(obj)
            % Run pipeline with interactive checkpoints
            
            status = 0;
            output = '';
            
            % Build initial command
            obj.command.clear();
            
            lastCheckpoint = 0;
            for i = 1:length(obj.checkpoints)
                checkpoint = obj.checkpoints{i};
                
                % Build command up to checkpoint
                for j = lastCheckpoint+1:checkpoint.index
                    obj.addTaskToCommand(j);
                end
                
                % Add write for checkpoint
                tempFile = obj.createTempFile(sprintf('checkpoint_%02d.h5', i));
                obj.command.addTask('write', 'filename', tempFile);
                
                % Execute
                fprintf('Executing to checkpoint: %s\n', checkpoint.name);
                [status, output] = obj.command.execute();
                
                if status ~= 0
                    error('Pipeline:ExecutionFailed', ...
                        'Execution failed at checkpoint %s', checkpoint.name);
                end
                
                % Load and visualize
                data = obj.loadDataFromFile(tempFile);
                
                if obj.enableCache
                    obj.cache.save(checkpoint.name, data);
                end
                
                % Visualize
                fprintf('Checkpoint reached: %s\n', checkpoint.name);
                obj.visualizeData(data, struct('title', checkpoint.name));
                
                % Custom visualization if provided
                if ~isempty(checkpoint.visualize)
                    checkpoint.visualize(data);
                end
                
                notify(obj, 'CheckpointReached');
                
                % Wait for user
                response = questdlg('Continue pipeline?', 'Checkpoint', ...
                    'Continue', 'Stop', 'Save & Stop', 'Continue');
                
                if strcmp(response, 'Stop')
                    return;
                elseif strcmp(response, 'Save & Stop')
                    obj.saveState();
                    return;
                end
                
                % Start new chain from checkpoint
                obj.command.clear();
                obj.command.addTask('read', 'path', tempFile);
                
                lastCheckpoint = checkpoint.index;
            end
            
            % Execute remaining tasks
            if lastCheckpoint < obj.length
                for j = lastCheckpoint+1:obj.length
                    obj.addTaskToCommand(j);
                end
                
                [status, output] = obj.command.execute();
            end
        end
        
        function addTaskToCommand(obj, index)
            % Add task at index to command
            
            task = obj.tasks{index};
            
            if iscell(task)
                % Parallel branch
                obj.command.addTaskChain(task);
            else
                % Single task
                obj.command.addTask(task.name, ...
                    obj.structToVarargin(task.getParameters()));
            end
        end
        
        function filename = createTempFile(obj, name)
            % Create temporary file
            
            tempDir = obj.config.get('output.checkpoints', tempdir);
            filename = fullfile(tempDir, name);
            obj.tempFiles{end+1} = filename;
        end
        
        function data = loadDataFromFile(obj, filename)
            % Load data from file
            
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case {'.h5', '.hdf5'}
                    % Find main dataset
                    info = h5info(filename);
                    datasets = obj.findDatasets(info);
                    if ~isempty(datasets)
                        data = h5read(filename, datasets{1});
                    else
                        error('Pipeline:NoData', 'No datasets found in file');
                    end
                    
                case {'.tif', '.tiff'}
                    data = imread(filename);
                    
                otherwise
                    error('Pipeline:UnsupportedFormat', ...
                        'Unsupported file format: %s', ext);
            end
        end
        
        function datasets = findDatasets(obj, info, prefix)
            % Find datasets in HDF5 info structure
            
            if nargin < 3
                prefix = '';
            end
            
            datasets = {};
            
            % Check datasets at current level
            for i = 1:length(info.Datasets)
                datasetPath = fullfile(prefix, info.Datasets(i).Name);
                datasets{end+1} = datasetPath;
            end
            
            % Check groups
            for i = 1:length(info.Groups)
                groupDatasets = obj.findDatasets(info.Groups(i), ...
                    info.Groups(i).Name);
                datasets = [datasets, groupDatasets];
            end
        end
        
        function tf = isSinogram(obj, data)
            % Check if data is sinogram
            tf = false;
            
            % Simple heuristic - customize as needed
            if ndims(data) == 2
                [h, w] = size(data);
                if h < w * 0.5  % More projections than detector pixels
                    tf = true;
                end
            end
        end
        
        function tf = isProjection(obj, data)
            % Check if data is projection
            tf = ndims(data) == 2 || ndims(data) == 3;
        end
        
        function tf = isVolume(obj, data)
            % Check if data is volume
            tf = ndims(data) == 3;
        end
        
        function saveState(obj)
            % Save pipeline state
            
            state = struct();
            state.config = obj.config;
            state.checkpoints = obj.checkpoints;
            state.tasks = obj.tasks;
            
            stateFile = fullfile(obj.config.get('output.checkpoints'), ...
                sprintf('pipeline_state_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
            
            save(stateFile, 'state');
            fprintf('Pipeline state saved to: %s\n', stateFile);
        end
        
        function cleanupTempFiles(obj)
            % Clean up temporary files
            
            for i = 1:length(obj.tempFiles)
                if exist(obj.tempFiles{i}, 'file')
                    delete(obj.tempFiles{i});
                end
            end
            
            obj.tempFiles = {};
        end
        
        function saveProfilingReport(obj, profData)
            % Save profiling report
            
            reportFile = fullfile(obj.config.get('output.logs'), ...
                sprintf('profile_report_%s.html', datestr(now, 'yyyymmdd_HHMMSS')));
            
            profsave(profData, reportFile);
            fprintf('Profiling report saved to: %s\n', reportFile);
        end
        
        function args = structToVarargin(obj, s)
            % Convert struct to varargin
            fields = fieldnames(s);
            args = cell(1, 2 * length(fields));
            
            for i = 1:length(fields)
                args{2*i-1} = fields{i};
                args{2*i} = s.(fields{i});
            end
        end
    end
end