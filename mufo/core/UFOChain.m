classdef UFOChain < handle
    % UFOChain - Chain multiple UFO tasks together
    %
    % This class provides functionality to build and execute chains of
    % UFO tasks with support for parallel branches and validation.
    %
    % Example:
    %   chain = mufo.core.UFOChain();
    %   chain.addTask('read', 'path', 'data.h5');
    %   chain.addTask('transpose');
    %   chain.addParallel({...
    %       mufo.core.UFOTask.create('fft', 'dimensions', 1), ...
    %       mufo.core.UFOTask.create('filter', 'type', 'ramp')});
    %   chain.execute();
    
    properties (Access = private)
        tasks           % Cell array of tasks/chains
        command         % UFOCommand instance
        cache           % Data cache for checkpoints
        validator       % Chain validator
    end
    
    properties (SetAccess = private)
        length          % Number of tasks in chain
        isValid         % Whether chain is valid
        executable      % UFO executable to use
    end
    
    properties
        name            % Chain name (optional)
        description     % Chain description (optional)
        saveCheckpoints % Save intermediate results
        checkpointDir   % Directory for checkpoints
        verbose         % Verbose output
    end
    
    events
        TaskAdded
        ChainExecuted
        ValidationFailed
    end
    
    methods
        function obj = UFOChain(varargin)
            % Constructor
            %
            % Parameters:
            %   'Name' - Chain name
            %   'Description' - Chain description
            %   'Executable' - UFO executable (default: 'ufo-launch')
            %   'SaveCheckpoints' - Save intermediate results (default: false)
            %   'CheckpointDir' - Directory for checkpoints
            %   'Verbose' - Verbose output (default: false)
            
            p = inputParser;
            addParameter(p, 'Name', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Description', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Executable', 'ufo-launch', @(x) ischar(x) || isstring(x));
            addParameter(p, 'SaveCheckpoints', false, @islogical);
            addParameter(p, 'CheckpointDir', tempdir, @(x) ischar(x) || isstring(x));
            addParameter(p, 'Verbose', false, @islogical);
            parse(p, varargin{:});
            
            obj.name = p.Results.Name;
            obj.description = p.Results.Description;
            obj.executable = p.Results.Executable;
            obj.saveCheckpoints = p.Results.SaveCheckpoints;
            obj.checkpointDir = p.Results.CheckpointDir;
            obj.verbose = p.Results.Verbose;
            
            obj.tasks = {};
            obj.length = 0;
            
            % Create command builder
            obj.command = mufo.core.UFOCommand(obj.executable, ...
                'Verbose', obj.verbose);
        end
        
        function addTask(obj, task, varargin)
            % Add a task to the chain
            %
            % Inputs:
            %   task - UFOTask object, task name, or task struct
            %   varargin - Parameters if task is a name
            
            if isa(task, 'mufo.core.UFOTask')
                % Direct task object
                taskObj = task;
            elseif ischar(task) || isstring(task)
                % Task name with parameters
                taskObj = mufo.core.UFOTask.create(task, varargin{:});
            elseif isstruct(task)
                % Task struct
                taskObj = obj.structToTask(task);
            else
                error('UFOChain:InvalidTask', ...
                    'Task must be UFOTask object, task name, or struct');
            end
            
            % Validate task
            if ~taskObj.isValid
                missing = taskObj.getMissingParameters();
                error('UFOChain:InvalidTask', ...
                    'Task "%s" is missing required parameters: %s', ...
                    taskObj.name, strjoin(missing, ', '));
            end
            
            obj.tasks{end+1} = taskObj;
            obj.length = obj.length + 1;
            
            notify(obj, 'TaskAdded');
        end
        
        function addParallel(obj, tasks)
            % Add parallel task branch
            %
            % Input:
            %   tasks - Cell array of tasks to run in parallel
            
            if ~iscell(tasks) || isempty(tasks)
                error('UFOChain:InvalidInput', ...
                    'Parallel tasks must be a non-empty cell array');
            end
            
            % Convert all tasks to UFOTask objects
            parallelTasks = cell(size(tasks));
            for i = 1:length(tasks)
                if isa(tasks{i}, 'mufo.core.UFOTask')
                    parallelTasks{i} = tasks{i};
                elseif ischar(tasks{i}) || isstring(tasks{i})
                    parallelTasks{i} = mufo.core.UFOTask.create(tasks{i});
                elseif isstruct(tasks{i})
                    parallelTasks{i} = obj.structToTask(tasks{i});
                else
                    error('UFOChain:InvalidTask', ...
                        'Invalid task type in parallel branch');
                end
                
                % Validate
                if ~parallelTasks{i}.isValid
                    error('UFOChain:InvalidTask', ...
                        'Task %d in parallel branch is not valid', i);
                end
            end
            
            obj.tasks{end+1} = parallelTasks;
            obj.length = obj.length + 1;
            
            notify(obj, 'TaskAdded');
        end
        
        function insertTask(obj, index, task, varargin)
            % Insert a task at specific position
            %
            % Inputs:
            %   index - Position to insert (1-based)
            %   task - Task to insert
            %   varargin - Task parameters if task is a name
            
            if index < 1 || index > obj.length + 1
                error('UFOChain:InvalidIndex', ...
                    'Index out of bounds');
            end
            
            % Create task object
            if isa(task, 'mufo.core.UFOTask')
                taskObj = task;
            elseif ischar(task) || isstring(task)
                taskObj = mufo.core.UFOTask.create(task, varargin{:});
            else
                error('UFOChain:InvalidTask', ...
                    'Task must be UFOTask object or task name');
            end
            
            % Insert into chain
            if index == 1
                obj.tasks = [{taskObj}, obj.tasks];
            elseif index > obj.length
                obj.tasks{end+1} = taskObj;
            else
                obj.tasks = [obj.tasks(1:index-1), {taskObj}, obj.tasks(index:end)];
            end
            
            obj.length = obj.length + 1;
        end
        
        function removeTask(obj, index)
            % Remove a task from the chain
            %
            % Input:
            %   index - Task index to remove
            
            if index < 1 || index > obj.length
                error('UFOChain:InvalidIndex', ...
                    'Index out of bounds');
            end
            
            obj.tasks(index) = [];
            obj.length = obj.length - 1;
        end
        
        function task = getTask(obj, index)
            % Get a task from the chain
            %
            % Input:
            %   index - Task index
            % Output:
            %   task - Task object or parallel task array
            
            if index < 1 || index > obj.length
                error('UFOChain:InvalidIndex', ...
                    'Index out of bounds');
            end
            
            task = obj.tasks{index};
        end
        
        function [status, output] = execute(obj, varargin)
            % Execute the task chain
            %
            % Parameters:
            %   'Async' - Execute asynchronously (default: false)
            %   'OutputFile' - Final output file (optional)
            %   'Process' - Use UFOProcess for monitoring (default: false)
            % Outputs:
            %   status - Exit status
            %   output - Command output
            
            p = inputParser;
            addParameter(p, 'Async', false, @islogical);
            addParameter(p, 'OutputFile', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'Process', false, @islogical);
            parse(p, varargin{:});
            
            % Validate chain
            if ~obj.validate()
                error('UFOChain:InvalidChain', ...
                    'Chain validation failed');
            end
            
            % Build command
            obj.buildCommand();
            
            % Add final output if specified
            if ~isempty(p.Results.OutputFile)
                obj.command.addTask('write', 'filename', p.Results.OutputFile);
            end
            
            % Execute
            if p.Results.Process
                % Use UFOProcess for monitoring
                proc = mufo.core.UFOProcess();
                proc.setCommand(obj.command);
                
                if p.Results.Async
                    proc.executeAsync();
                    status = 0;
                    output = proc;  % Return process handle
                else
                    [status, output] = proc.execute();
                end
            else
                % Direct execution
                [status, output] = obj.command.execute();
            end
            
            notify(obj, 'ChainExecuted');
        end
        
        function cmd = buildCommandString(obj)
            % Build and return the command string
            %
            % Output:
            %   cmd - Complete UFO command string
            
            obj.buildCommand();
            cmd = obj.command.buildCommand();
        end
        
        function valid = validate(obj)
            % Validate the chain
            %
            % Output:
            %   valid - Whether chain is valid
            
            valid = true;
            
            if obj.length == 0
                warning('UFOChain:EmptyChain', 'Chain is empty');
                valid = false;
                return;
            end
            
            % Validate each task
            for i = 1:obj.length
                if iscell(obj.tasks{i})
                    % Parallel branch
                    for j = 1:length(obj.tasks{i})
                        if ~obj.tasks{i}{j}.isValid
                            valid = false;
                            notify(obj, 'ValidationFailed');
                            return;
                        end
                    end
                else
                    % Single task
                    if ~obj.tasks{i}.isValid
                        valid = false;
                        notify(obj, 'ValidationFailed');
                        return;
                    end
                end
            end
            
            % Additional validation can be added here
            % e.g., checking data flow compatibility
            
            obj.isValid = valid;
        end
        
        function clear(obj)
            % Clear all tasks from the chain
            obj.tasks = {};
            obj.length = 0;
            obj.command.clear();
        end
        
        function info = getChainInfo(obj)
            % Get information about the chain
            info = struct();
            info.name = obj.name;
            info.description = obj.description;
            info.length = obj.length;
            info.isValid = obj.validate();
            info.tasks = cell(obj.length, 1);
            
            for i = 1:obj.length
                if iscell(obj.tasks{i})
                    % Parallel branch
                    info.tasks{i} = struct('type', 'parallel', ...
                        'count', length(obj.tasks{i}));
                else
                    % Single task
                    info.tasks{i} = obj.tasks{i}.getTaskInfo();
                end
            end
        end
        
        function saveChain(obj, filename)
            % Save chain configuration to file
            %
            % Input:
            %   filename - Output filename (.mat or .json)
            
            [~, ~, ext] = fileparts(filename);
            
            % Create chain structure
            chainStruct = struct();
            chainStruct.name = obj.name;
            chainStruct.description = obj.description;
            chainStruct.executable = obj.executable;
            chainStruct.tasks = obj.tasksToStruct();
            
            switch lower(ext)
                case '.mat'
                    save(filename, 'chainStruct');
                case '.json'
                    jsonStr = jsonencode(chainStruct, 'PrettyPrint', true);
                    fid = fopen(filename, 'w');
                    fprintf(fid, '%s', jsonStr);
                    fclose(fid);
                otherwise
                    error('UFOChain:UnsupportedFormat', ...
                        'Unsupported file format: %s', ext);
            end
        end
        
        function loadChain(obj, filename)
            % Load chain configuration from file
            %
            % Input:
            %   filename - Input filename (.mat or .json)
            
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case '.mat'
                    data = load(filename);
                    chainStruct = data.chainStruct;
                case '.json'
                    jsonStr = fileread(filename);
                    chainStruct = jsondecode(jsonStr);
                otherwise
                    error('UFOChain:UnsupportedFormat', ...
                        'Unsupported file format: %s', ext);
            end
            
            % Clear current chain
            obj.clear();
            
            % Load configuration
            obj.name = chainStruct.name;
            obj.description = chainStruct.description;
            
            % Rebuild tasks
            obj.structToTasks(chainStruct.tasks);
        end
    end
    
    methods (Access = private)
        function buildCommand(obj)
            % Build UFO command from task chain
            
            obj.command.clear();
            
            for i = 1:obj.length
                if iscell(obj.tasks{i})
                    % Parallel branch
                    parallelStructs = cell(size(obj.tasks{i}));
                    for j = 1:length(obj.tasks{i})
                        task = obj.tasks{i}{j};
                        taskStruct = struct('name', task.name, ...
                            'parameters', task.getParameters());
                        parallelStructs{j} = taskStruct;
                    end
                    obj.command.addTaskChain(parallelStructs);
                else
                    % Single task
                    task = obj.tasks{i};
                    obj.command.addTask(task.name, ...
                        obj.structToVarargin(task.getParameters()));
                end
                
                % Add checkpoint if needed
                if obj.saveCheckpoints && i < obj.length
                    checkpointFile = fullfile(obj.checkpointDir, ...
                        sprintf('checkpoint_%03d.h5', i));
                    obj.addCheckpoint(checkpointFile);
                end
            end
        end
        
        function addCheckpoint(obj, filename)
            % Add checkpoint write/read tasks
            obj.command.addTask('write', 'filename', filename);
            obj.command.addTask('read', 'path', filename);
        end
        
        function taskStruct = tasksToStruct(obj)
            % Convert tasks to structure array
            taskStruct = cell(obj.length, 1);
            
            for i = 1:obj.length
                if iscell(obj.tasks{i})
                    % Parallel branch
                    parallel = cell(length(obj.tasks{i}), 1);
                    for j = 1:length(obj.tasks{i})
                        task = obj.tasks{i}{j};
                        parallel{j} = struct('name', task.name, ...
                            'parameters', task.getParameters());
                    end
                    taskStruct{i} = struct('type', 'parallel', 'tasks', {parallel});
                else
                    % Single task
                    task = obj.tasks{i};
                    taskStruct{i} = struct('type', 'single', ...
                        'name', task.name, ...
                        'parameters', task.getParameters());
                end
            end
        end
        
        function structToTasks(obj, taskStruct)
            % Convert structure array to tasks
            for i = 1:length(taskStruct)
                if strcmp(taskStruct{i}.type, 'parallel')
                    % Parallel branch
                    parallel = cell(length(taskStruct{i}.tasks), 1);
                    for j = 1:length(taskStruct{i}.tasks)
                        taskData = taskStruct{i}.tasks{j};
                        parallel{j} = mufo.core.UFOTask.create(taskData.name);
                        
                        % Set parameters
                        params = fieldnames(taskData.parameters);
                        for k = 1:length(params)
                            parallel{j}.setParameter(params{k}, ...
                                taskData.parameters.(params{k}));
                        end
                    end
                    obj.addParallel(parallel);
                else
                    % Single task
                    task = mufo.core.UFOTask.create(taskStruct{i}.name);
                    
                    % Set parameters
                    params = fieldnames(taskStruct{i}.parameters);
                    for j = 1:length(params)
                        task.setParameter(params{j}, ...
                            taskStruct{i}.parameters.(params{j}));
                    end
                    
                    obj.addTask(task);
                end
            end
        end
        
        function task = structToTask(obj, taskStruct)
            % Convert struct to task object
            if ~isfield(taskStruct, 'name')
                error('UFOChain:InvalidStruct', ...
                    'Task struct must have a name field');
            end
            
            task = mufo.core.UFOTask.create(taskStruct.name);
            
            if isfield(taskStruct, 'parameters')
                params = fieldnames(taskStruct.parameters);
                for i = 1:length(params)
                    task.setParameter(params{i}, ...
                        taskStruct.parameters.(params{i}));
                end
            end
        end
        
        function args = structToVarargin(obj, s)
            % Convert struct to varargin cell array
            fields = fieldnames(s);
            args = cell(1, 2 * length(fields));
            
            for i = 1:length(fields)
                args{2*i-1} = fields{i};
                args{2*i} = s.(fields{i});
            end
        end
    end
end