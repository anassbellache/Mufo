classdef UFOCommand < handle
    % UFOCommand - Base class for building UFO command-line operations
    %
    % This class provides the foundation for constructing UFO commands
    % with proper parameter handling, validation, and execution.
    %
    % Example:
    %   cmd = mufo.core.UFOCommand('ufo-launch');
    %   cmd.addTask('read', 'path', 'data.h5');
    %   cmd.addTask('write', 'filename', 'output.tif');
    %   cmd.execute();
    
    properties (Access = protected)
        executable      % UFO executable path
        tasks           % Cell array of task structures
        environment     % Environment variables
        workingDir      % Working directory for execution
        verbose         % Verbose output flag
        dryRun          % Dry run mode (build command without executing)
    end
    
    properties (Constant, Access = private)
        % Default UFO paths
        DEFAULT_UFO_BIN = '/usr/bin/ufo-launch';
        DEFAULT_PLUGIN_DIR = '/usr/lib/x86_64-linux-gnu/ufo';
        DEFAULT_KERNEL_PATH = '/usr/share/ufo/kernels';
    end
    
    methods
        function obj = UFOCommand(executable, varargin)
            % Constructor - Initialize UFO command builder
            %
            % Inputs:
            %   executable - UFO executable name or path (default: 'ufo-launch')
            %   'WorkingDir' - Working directory for execution
            %   'Verbose' - Enable verbose output (default: false)
            %   'DryRun' - Build command without executing (default: false)
            %   'Environment' - Struct with environment variables
            
            p = inputParser;
            addRequired(p, 'executable', @(x) ischar(x) || isstring(x));
            addParameter(p, 'WorkingDir', pwd, @(x) ischar(x) || isstring(x));
            addParameter(p, 'Verbose', false, @islogical);
            addParameter(p, 'DryRun', false, @islogical);
            addParameter(p, 'Environment', struct(), @isstruct);
            parse(p, executable, varargin{:});
            
            % Set executable path
            if strcmpi(executable, 'ufo-launch')
                obj.executable = obj.DEFAULT_UFO_BIN;
            elseif strcmpi(executable, 'tofu')
                obj.executable = 'tofu';
            else
                obj.executable = executable;
            end
            
            % Validate executable exists
            if ~obj.isExecutableAvailable()
                warning('UFOCommand:ExecutableNotFound', ...
                    'UFO executable not found: %s', obj.executable);
            end
            
            obj.tasks = {};
            obj.workingDir = p.Results.WorkingDir;
            obj.verbose = p.Results.Verbose;
            obj.dryRun = p.Results.DryRun;
            
            % Setup environment
            obj.environment = obj.getDefaultEnvironment();
            envFields = fieldnames(p.Results.Environment);
            for i = 1:length(envFields)
                obj.environment.(envFields{i}) = p.Results.Environment.(envFields{i});
            end
        end
        
        function addTask(obj, taskName, varargin)
            % Add a task to the command chain
            %
            % Inputs:
            %   taskName - Name of the UFO task
            %   varargin - Parameter name-value pairs
            
            task = struct();
            task.name = taskName;
            task.parameters = struct();
            
            % Parse parameters
            if mod(length(varargin), 2) ~= 0
                error('UFOCommand:InvalidParameters', ...
                    'Parameters must be in name-value pairs');
            end
            
            for i = 1:2:length(varargin)
                paramName = varargin{i};
                paramValue = varargin{i+1};
                
                % Validate parameter name
                if ~ischar(paramName) && ~isstring(paramName)
                    error('UFOCommand:InvalidParameterName', ...
                        'Parameter names must be strings');
                end
                
                % Convert parameter value to string
                task.parameters.(paramName) = obj.convertToString(paramValue);
            end
            
            obj.tasks{end+1} = task;
        end
        
        function addTaskChain(obj, tasks)
            % Add multiple tasks as a parallel chain
            %
            % Inputs:
            %   tasks - Cell array of task structures or UFOTask objects
            
            if ~iscell(tasks)
                error('UFOCommand:InvalidTaskChain', ...
                    'Task chain must be a cell array');
            end
            
            obj.tasks{end+1} = tasks;
        end
        
        function cmdStr = buildCommand(obj)
            % Build the complete command string
            
            if isempty(obj.tasks)
                error('UFOCommand:NoTasks', ...
                    'No tasks added to command');
            end
            
            cmdStr = obj.executable;
            
            % Add tasks
            for i = 1:length(obj.tasks)
                if i > 1
                    cmdStr = [cmdStr, ' ! '];
                else
                    cmdStr = [cmdStr, ' '];
                end
                
                if iscell(obj.tasks{i})
                    % Parallel task chain
                    cmdStr = [cmdStr, obj.buildParallelChain(obj.tasks{i})];
                else
                    % Single task
                    cmdStr = [cmdStr, obj.buildTaskString(obj.tasks{i})];
                end
            end
        end
        
        function [status, output] = execute(obj)
            % Execute the UFO command
            %
            % Outputs:
            %   status - Exit status (0 = success)
            %   output - Command output string
            
            % Build command
            cmdStr = obj.buildCommand();
            
            if obj.verbose || obj.dryRun
                fprintf('UFO Command: %s\n', cmdStr);
            end
            
            if obj.dryRun
                status = 0;
                output = sprintf('Dry run - command not executed:\n%s', cmdStr);
                return;
            end
            
            % Setup environment
            obj.setupEnvironment();
            
            % Change to working directory
            oldDir = pwd;
            cd(obj.workingDir);
            
            try
                % Execute command
                if obj.verbose
                    [status, output] = system(cmdStr, '-echo');
                else
                    [status, output] = system(cmdStr);
                end
                
                if status ~= 0
                    warning('UFOCommand:ExecutionFailed', ...
                        'UFO command failed with status %d:\n%s', status, output);
                end
            catch ME
                cd(oldDir);
                rethrow(ME);
            end
            
            cd(oldDir);
        end
        
        function clear(obj)
            % Clear all tasks
            obj.tasks = {};
        end
        
        function setEnvironment(obj, name, value)
            % Set an environment variable
            obj.environment.(name) = value;
        end
        
        function env = getEnvironment(obj)
            % Get current environment settings
            env = obj.environment;
        end
    end
    
    methods (Access = protected)
        function env = getDefaultEnvironment(obj)
            % Get default UFO environment variables
            env = struct();
            env.UFO_PLUGIN_DIR = obj.DEFAULT_PLUGIN_DIR;
            env.UFO_KERNEL_PATH = obj.DEFAULT_KERNEL_PATH;
            env.G_MESSAGES_DEBUG = '';
            env.UFO_DEVICE_TYPE = 'gpu';
            env.UFO_DISABLE_OPENCL_CACHE = '0';
        end
        
        function setupEnvironment(obj)
            % Set environment variables for UFO execution
            fields = fieldnames(obj.environment);
            for i = 1:length(fields)
                value = obj.environment.(fields{i});
                if ~isempty(value)
                    setenv(fields{i}, value);
                end
            end
        end
        
        function str = buildTaskString(obj, task)
            % Build string representation of a single task
            str = task.name;
            
            % Add parameters
            params = fieldnames(task.parameters);
            for i = 1:length(params)
                str = sprintf('%s %s=%s', str, params{i}, task.parameters.(params{i}));
            end
        end
        
        function str = buildParallelChain(obj, tasks)
            % Build string for parallel task chain
            str = '[ ';
            
            for i = 1:length(tasks)
                if i > 1
                    str = [str, ', '];
                end
                
                if isstruct(tasks{i})
                    str = [str, obj.buildTaskString(tasks{i})];
                elseif isa(tasks{i}, 'mufo.core.UFOTask')
                    str = [str, tasks{i}.buildCommand()];
                else
                    error('UFOCommand:InvalidTask', ...
                        'Invalid task type in parallel chain');
                end
            end
            
            str = [str, ' ]'];
        end
        
        function str = convertToString(obj, value)
            % Convert parameter value to string representation
            if ischar(value) || isstring(value)
                str = char(value);
            elseif isnumeric(value)
                if isscalar(value)
                    if mod(value, 1) == 0
                        str = sprintf('%d', value);
                    else
                        str = sprintf('%.6g', value);
                    end
                else
                    % Array values
                    str = sprintf('%.6g,', value);
                    str(end) = []; % Remove trailing comma
                end
            elseif islogical(value)
                if value
                    str = 'TRUE';
                else
                    str = 'FALSE';
                end
            else
                error('UFOCommand:UnsupportedType', ...
                    'Unsupported parameter type: %s', class(value));
            end
        end
        
        function available = isExecutableAvailable(obj)
            % Check if UFO executable is available
            if ispc
                cmd = sprintf('where %s', obj.executable);
            else
                cmd = sprintf('which %s', obj.executable);
            end
            
            [status, ~] = system(cmd);
            available = (status == 0);
        end
    end
end