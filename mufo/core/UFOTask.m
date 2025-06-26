classdef UFOTask < handle
    % UFOTask - Base class for UFO task wrappers
    %
    % This class provides a foundation for wrapping individual UFO tasks
    % with parameter validation, type checking, and command generation.
    %
    % Example:
    %   task = mufo.core.UFOTask('read');
    %   task.setParameter('path', 'data.h5');
    %   task.setParameter('y', 100);
    %   cmd = task.buildCommand();
    
    properties (Access = protected)
        name            % Task name
        parameters      % Task parameters
        parameterDefs   % Parameter definitions
        required        % Required parameters
        validated       % Validation flag
    end
    
    properties (Dependent)
        isValid         % Whether task configuration is valid
    end
    
    methods
        function obj = UFOTask(name)
            % Constructor
            %
            % Input:
            %   name - UFO task name
            
            if nargin < 1
                error('UFOTask:NoName', 'Task name is required');
            end
            
            obj.name = name;
            obj.parameters = struct();
            obj.parameterDefs = struct();
            obj.required = {};
            obj.validated = false;
            
            % Define parameters for this task
            obj.defineParameters();
        end
        
        function setParameter(obj, param, value)
            % Set a task parameter
            %
            % Inputs:
            %   param - Parameter name
            %   value - Parameter value
            
            % Check if parameter is defined
            if ~isfield(obj.parameterDefs, param)
                warning('UFOTask:UnknownParameter', ...
                    'Unknown parameter "%s" for task "%s"', param, obj.name);
            else
                % Validate parameter value
                if obj.validateParameter(param, value)
                    obj.parameters.(param) = value;
                    obj.validated = false; % Need revalidation
                else
                    error('UFOTask:InvalidValue', ...
                        'Invalid value for parameter "%s"', param);
                end
            end
        end
        
        function value = getParameter(obj, param)
            % Get a parameter value
            %
            % Input:
            %   param - Parameter name
            % Output:
            %   value - Parameter value (empty if not set)
            
            if isfield(obj.parameters, param)
                value = obj.parameters.(param);
            else
                value = [];
            end
        end
        
        function params = getParameters(obj)
            % Get all parameters
            params = obj.parameters;
        end
        
        function setParameters(obj, varargin)
            % Set multiple parameters at once
            %
            % Input:
            %   varargin - Name-value pairs
            
            if mod(length(varargin), 2) ~= 0
                error('UFOTask:InvalidInput', ...
                    'Parameters must be in name-value pairs');
            end
            
            for i = 1:2:length(varargin)
                obj.setParameter(varargin{i}, varargin{i+1});
            end
        end
        
        function cmd = buildCommand(obj)
            % Build command string for this task
            
            % Validate task first
            if ~obj.isValid
                error('UFOTask:InvalidConfiguration', ...
                    'Task configuration is not valid. Missing required parameters.');
            end
            
            cmd = obj.name;
            
            % Add parameters
            params = fieldnames(obj.parameters);
            for i = 1:length(params)
                value = obj.parameters.(params{i});
                valueStr = obj.convertToString(value);
                cmd = sprintf('%s %s=%s', cmd, params{i}, valueStr);
            end
        end
        
        function valid = get.isValid(obj)
            % Check if task configuration is valid
            
            % Check all required parameters are set
            for i = 1:length(obj.required)
                if ~isfield(obj.parameters, obj.required{i})
                    valid = false;
                    return;
                end
            end
            
            % Run custom validation if implemented
            valid = obj.validate();
            obj.validated = valid;
        end
        
        function missing = getMissingParameters(obj)
            % Get list of missing required parameters
            missing = {};
            
            for i = 1:length(obj.required)
                if ~isfield(obj.parameters, obj.required{i})
                    missing{end+1} = obj.required{i};
                end
            end
        end
        
        function info = getTaskInfo(obj)
            % Get task information
            info = struct();
            info.name = obj.name;
            info.parameters = obj.parameters;
            info.required = obj.required;
            info.isValid = obj.isValid;
            info.missing = obj.getMissingParameters();
            
            % Add parameter definitions
            info.parameterDefs = obj.parameterDefs;
        end
        
        function reset(obj)
            % Reset all parameters
            obj.parameters = struct();
            obj.validated = false;
        end
        
        function clone = copy(obj)
            % Create a copy of this task
            clone = feval(class(obj), obj.name);
            
            % Copy parameters
            params = fieldnames(obj.parameters);
            for i = 1:length(params)
                clone.parameters.(params{i}) = obj.parameters.(params{i});
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define parameters for this task
            % Override in subclasses for specific tasks
            
            % This is a generic implementation
            % Specific task classes should override this
        end
        
        function valid = validate(obj)
            % Custom validation logic
            % Override in subclasses for specific validation
            valid = true;
        end
        
        function valid = validateParameter(obj, param, value)
            % Validate a parameter value
            %
            % Inputs:
            %   param - Parameter name
            %   value - Parameter value
            % Output:
            %   valid - Whether value is valid
            
            valid = true;
            
            % Check against parameter definition
            if isfield(obj.parameterDefs, param)
                def = obj.parameterDefs.(param);
                
                % Type checking
                if isfield(def, 'type')
                    switch def.type
                        case 'numeric'
                            valid = isnumeric(value);
                        case 'integer'
                            valid = isnumeric(value) && all(mod(value, 1) == 0);
                        case 'positive'
                            valid = isnumeric(value) && all(value > 0);
                        case 'string'
                            valid = ischar(value) || isstring(value);
                        case 'logical'
                            valid = islogical(value);
                        case 'file'
                            valid = (ischar(value) || isstring(value));
                        case 'directory'
                            valid = (ischar(value) || isstring(value));
                    end
                end
                
                % Range checking
                if valid && isfield(def, 'range') && isnumeric(value)
                    valid = all(value >= def.range(1) & value <= def.range(2));
                end
                
                % Enumeration checking
                if valid && isfield(def, 'values')
                    if ischar(value) || isstring(value)
                        valid = any(strcmpi(value, def.values));
                    else
                        valid = any(value == def.values);
                    end
                end
            end
        end
        
        function str = convertToString(obj, value)
            % Convert parameter value to string
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
                error('UFOTask:UnsupportedType', ...
                    'Cannot convert %s to string', class(value));
            end
        end
        
        function addParameterDef(obj, name, type, varargin)
            % Add a parameter definition
            %
            % Inputs:
            %   name - Parameter name
            %   type - Parameter type
            %   varargin - Additional properties (range, values, etc.)
            
            p = inputParser;
            addRequired(p, 'name', @(x) ischar(x) || isstring(x));
            addRequired(p, 'type', @(x) ischar(x) || isstring(x));
            addParameter(p, 'default', []);
            addParameter(p, 'range', []);
            addParameter(p, 'values', {});
            addParameter(p, 'description', '');
            addParameter(p, 'required', false, @islogical);
            parse(p, name, type, varargin{:});
            
            def = struct();
            def.type = p.Results.type;
            def.default = p.Results.default;
            def.range = p.Results.range;
            def.values = p.Results.values;
            def.description = p.Results.description;
            
            obj.parameterDefs.(name) = def;
            
            % Add to required list if needed
            if p.Results.required
                obj.required{end+1} = name;
            end
            
            % Set default value if provided
            if ~isempty(def.default)
                obj.parameters.(name) = def.default;
            end
        end
    end
    
    methods (Static)
        function task = create(taskName, varargin)
            % Factory method to create specific task instances
            %
            % Input:
            %   taskName - Name of the task
            %   varargin - Parameters to set
            % Output:
            %   task - Task instance
            
            % Try to find specific task class
            className = sprintf('mufo.tasks.%s', ...
                strcat(upper(taskName(1)), taskName(2:end)));
            
            if exist(className, 'class')
                % Use specific task class
                task = feval(className);
            else
                % Use generic task
                task = mufo.core.UFOTask(taskName);
            end
            
            % Set parameters
            if ~isempty(varargin)
                task.setParameters(varargin{:});
            end
        end
    end
end