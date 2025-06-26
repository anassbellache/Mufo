classdef UFOProcess < handle
    % UFOProcess - UFO process executor with monitoring and control
    %
    % This class manages UFO process execution including asynchronous
    % operations, progress monitoring, and resource management.
    %
    % Example:
    %   proc = mufo.core.UFOProcess();
    %   proc.setCommand('ufo-launch read path=data.h5 ! write filename=out.tif');
    %   proc.executeAsync();
    %   while proc.isRunning()
    %       fprintf('Progress: %.1f%%\n', proc.getProgress());
    %       pause(1);
    %   end
    
    properties (Access = private)
        command         % Command to execute
        processHandle   % Handle to running process
        startTime       % Process start time
        endTime         % Process end time
        status          % Exit status
        output          % Process output
        errorOutput     % Error output
        logFile         % Log file path
        progressFile    % Progress monitoring file
        isAsync         % Asynchronous execution flag
        gpuMonitor      % GPU monitor handle
    end
    
    properties (SetAccess = private)
        state           % Current process state
        progress        % Current progress (0-100)
        gpuMemoryUsage  % GPU memory usage
        gpuUtilization  % GPU utilization percentage
    end
    
    properties (Constant)
        % Process states
        STATE_IDLE = 'idle';
        STATE_RUNNING = 'running';
        STATE_COMPLETED = 'completed';
        STATE_FAILED = 'failed';
        STATE_CANCELLED = 'cancelled';
    end
    
    events
        ProcessStarted
        ProcessCompleted
        ProcessFailed
        ProgressUpdate
    end
    
    methods
        function obj = UFOProcess(varargin)
            % Constructor
            %
            % Parameters:
            %   'LogFile' - Path to log file (optional)
            %   'MonitorGPU' - Enable GPU monitoring (default: true)
            
            p = inputParser;
            addParameter(p, 'LogFile', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'MonitorGPU', true, @islogical);
            parse(p, varargin{:});
            
            obj.state = obj.STATE_IDLE;
            obj.progress = 0;
            obj.logFile = p.Results.LogFile;
            obj.gpuMemoryUsage = 0;
            obj.gpuUtilization = 0;
            
            % Setup GPU monitoring if requested
            if p.Results.MonitorGPU && obj.isGPUAvailable()
                obj.setupGPUMonitor();
            end
        end
        
        function setCommand(obj, command)
            % Set the command to execute
            %
            % Input:
            %   command - UFO command string or UFOCommand object
            
            if isa(command, 'mufo.core.UFOCommand')
                obj.command = command.buildCommand();
            elseif ischar(command) || isstring(command)
                obj.command = char(command);
            else
                error('UFOProcess:InvalidCommand', ...
                    'Command must be a string or UFOCommand object');
            end
        end
        
        function [status, output] = execute(obj)
            % Execute command synchronously
            %
            % Outputs:
            %   status - Exit status (0 = success)
            %   output - Combined output string
            
            obj.validateReady();
            obj.isAsync = false;
            
            obj.startExecution();
            
            try
                % Execute command
                if ~isempty(obj.logFile)
                    logCmd = sprintf('%s 2>&1 | tee %s', obj.command, obj.logFile);
                    [obj.status, obj.output] = system(logCmd);
                else
                    [obj.status, obj.output] = system(obj.command);
                end
                
                obj.endTime = datetime('now');
                
                if obj.status == 0
                    obj.state = obj.STATE_COMPLETED;
                    notify(obj, 'ProcessCompleted');
                else
                    obj.state = obj.STATE_FAILED;
                    notify(obj, 'ProcessFailed');
                end
                
            catch ME
                obj.state = obj.STATE_FAILED;
                obj.endTime = datetime('now');
                notify(obj, 'ProcessFailed');
                rethrow(ME);
            end
            
            status = obj.status;
            output = obj.output;
        end
        
        function executeAsync(obj)
            % Execute command asynchronously
            
            obj.validateReady();
            obj.isAsync = true;
            
            obj.startExecution();
            
            % Create temporary files for output capture
            obj.setupTempFiles();
            
            % Build command with output redirection
            asyncCmd = obj.buildAsyncCommand();
            
            % Start process in background
            if isunix
                % Unix/Linux
                obj.processHandle = system(sprintf('%s &', asyncCmd));
                
                % Get process ID
                [~, pidStr] = system('echo $!');
                obj.processHandle = str2double(strtrim(pidStr));
            else
                % Windows
                error('UFOProcess:NotImplemented', ...
                    'Asynchronous execution not yet implemented for Windows');
            end
            
            % Start monitoring timer
            obj.startMonitoring();
        end
        
        function cancel(obj)
            % Cancel running process
            
            if ~strcmp(obj.state, obj.STATE_RUNNING)
                warning('UFOProcess:NotRunning', ...
                    'No process to cancel');
                return;
            end
            
            if isunix && ~isempty(obj.processHandle)
                % Send SIGTERM
                system(sprintf('kill -15 %d', obj.processHandle));
                pause(2);
                
                % Check if still running and force kill
                [status, ~] = system(sprintf('ps -p %d', obj.processHandle));
                if status == 0
                    system(sprintf('kill -9 %d', obj.processHandle));
                end
            end
            
            obj.state = obj.STATE_CANCELLED;
            obj.endTime = datetime('now');
        end
        
        function running = isRunning(obj)
            % Check if process is currently running
            running = strcmp(obj.state, obj.STATE_RUNNING);
            
            % Double-check for async processes
            if running && obj.isAsync && ~isempty(obj.processHandle)
                if isunix
                    [status, ~] = system(sprintf('ps -p %d', obj.processHandle));
                    if status ~= 0
                        % Process no longer exists
                        obj.checkCompletion();
                        running = false;
                    end
                end
            end
        end
        
        function wait(obj, timeout)
            % Wait for process completion
            %
            % Input:
            %   timeout - Maximum wait time in seconds (optional)
            
            if nargin < 2
                timeout = inf;
            end
            
            startWait = tic;
            
            while obj.isRunning()
                pause(0.1);
                
                if toc(startWait) > timeout
                    warning('UFOProcess:Timeout', ...
                        'Process wait timeout after %.1f seconds', timeout);
                    break;
                end
            end
        end
        
        function prog = getProgress(obj)
            % Get current progress (0-100)
            prog = obj.progress;
            
            % Try to update progress if running
            if obj.isRunning() && ~isempty(obj.progressFile)
                obj.updateProgress();
            end
        end
        
        function [memUsage, utilization] = getGPUStats(obj)
            % Get current GPU statistics
            memUsage = obj.gpuMemoryUsage;
            utilization = obj.gpuUtilization;
        end
        
        function duration = getExecutionTime(obj)
            % Get execution duration
            if isempty(obj.startTime)
                duration = 0;
            elseif isempty(obj.endTime)
                duration = seconds(datetime('now') - obj.startTime);
            else
                duration = seconds(obj.endTime - obj.startTime);
            end
        end
        
        function log = getOutput(obj)
            % Get process output
            log = obj.output;
            
            % Read from log file if available
            if ~isempty(obj.logFile) && exist(obj.logFile, 'file')
                log = fileread(obj.logFile);
            end
        end
        
        function info = getInfo(obj)
            % Get process information structure
            info = struct();
            info.command = obj.command;
            info.state = obj.state;
            info.progress = obj.progress;
            info.startTime = obj.startTime;
            info.endTime = obj.endTime;
            info.duration = obj.getExecutionTime();
            info.status = obj.status;
            info.gpuMemoryUsage = obj.gpuMemoryUsage;
            info.gpuUtilization = obj.gpuUtilization;
        end
    end
    
    methods (Access = private)
        function validateReady(obj)
            % Validate process is ready to execute
            if isempty(obj.command)
                error('UFOProcess:NoCommand', ...
                    'No command set for execution');
            end
            
            if obj.isRunning()
                error('UFOProcess:AlreadyRunning', ...
                    'Process is already running');
            end
        end
        
        function startExecution(obj)
            % Initialize execution
            obj.state = obj.STATE_RUNNING;
            obj.startTime = datetime('now');
            obj.endTime = [];
            obj.status = [];
            obj.output = '';
            obj.errorOutput = '';
            obj.progress = 0;
            
            notify(obj, 'ProcessStarted');
        end
        
        function setupTempFiles(obj)
            % Setup temporary files for async execution
            tempDir = tempname;
            mkdir(tempDir);
            
            obj.progressFile = fullfile(tempDir, 'progress.txt');
            
            if isempty(obj.logFile)
                obj.logFile = fullfile(tempDir, 'output.log');
            end
        end
        
        function cmd = buildAsyncCommand(obj)
            % Build command for asynchronous execution
            
            % Add progress monitoring if possible
            progressCmd = obj.addProgressMonitoring(obj.command);
            
            % Add output redirection
            if ~isempty(obj.logFile)
                cmd = sprintf('%s > %s 2>&1', progressCmd, obj.logFile);
            else
                cmd = progressCmd;
            end
        end
        
        function cmd = addProgressMonitoring(obj, baseCmd)
            % Add progress monitoring to command
            % This is a simplified version - real implementation would
            % parse UFO output for progress indicators
            
            cmd = baseCmd;
            
            % Try to detect if this is a reconstruction command
            if contains(baseCmd, 'backproject') || contains(baseCmd, 'ifft')
                % Wrap with progress monitoring
                cmd = sprintf('(%s) 2>&1 | tee %s', baseCmd, obj.progressFile);
            end
        end
        
        function startMonitoring(obj)
            % Start monitoring timer for async process
            obj.gpuMonitor = timer(...
                'Period', 1, ...
                'ExecutionMode', 'fixedRate', ...
                'TimerFcn', @(~,~) obj.updateMonitoring());
            
            start(obj.gpuMonitor);
        end
        
        function updateMonitoring(obj)
            % Update monitoring information
            if ~obj.isRunning()
                % Stop monitoring if process completed
                if ~isempty(obj.gpuMonitor) && isvalid(obj.gpuMonitor)
                    stop(obj.gpuMonitor);
                    delete(obj.gpuMonitor);
                end
                return;
            end
            
            % Update progress
            obj.updateProgress();
            
            % Update GPU stats
            obj.updateGPUStats();
            
            % Check if process completed
            if obj.isAsync && ~isempty(obj.processHandle)
                if isunix
                    [status, ~] = system(sprintf('ps -p %d', obj.processHandle));
                    if status ~= 0
                        obj.checkCompletion();
                    end
                end
            end
        end
        
        function updateProgress(obj)
            % Update progress from monitoring file
            if exist(obj.progressFile, 'file')
                try
                    % Parse progress from file
                    % This is simplified - real implementation would parse
                    % UFO-specific progress indicators
                    content = fileread(obj.progressFile);
                    
                    % Look for percentage indicators
                    matches = regexp(content, '(\d+)%', 'tokens');
                    if ~isempty(matches)
                        obj.progress = str2double(matches{end}{1});
                        notify(obj, 'ProgressUpdate');
                    end
                catch
                    % Ignore errors reading progress
                end
            end
        end
        
        function updateGPUStats(obj)
            % Update GPU statistics
            if obj.isGPUAvailable()
                try
                    % Use nvidia-smi to get GPU stats
                    [status, output] = system('nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader,nounits');
                    
                    if status == 0
                        values = sscanf(output, '%f, %f');
                        if length(values) >= 2
                            obj.gpuMemoryUsage = values(1); % MB
                            obj.gpuUtilization = values(2); % Percentage
                        end
                    end
                catch
                    % Ignore errors
                end
            end
        end
        
        function checkCompletion(obj)
            % Check and update completion status
            obj.endTime = datetime('now');
            
            % Read final output
            if exist(obj.logFile, 'file')
                obj.output = fileread(obj.logFile);
            end
            
            % Try to determine exit status
            % This is simplified - real implementation would check
            % for UFO-specific error patterns
            if contains(obj.output, 'error', 'IgnoreCase', true) || ...
               contains(obj.output, 'failed', 'IgnoreCase', true)
                obj.status = 1;
                obj.state = obj.STATE_FAILED;
                notify(obj, 'ProcessFailed');
            else
                obj.status = 0;
                obj.state = obj.STATE_COMPLETED;
                obj.progress = 100;
                notify(obj, 'ProcessCompleted');
            end
            
            % Cleanup monitoring
            if ~isempty(obj.gpuMonitor) && isvalid(obj.gpuMonitor)
                stop(obj.gpuMonitor);
                delete(obj.gpuMonitor);
            end
        end
        
        function setupGPUMonitor(obj)
            % Setup GPU monitoring capabilities
            % This is a placeholder - actual implementation would
            % integrate with GPU monitoring tools
        end
        
        function available = isGPUAvailable(obj)
            % Check if GPU monitoring is available
            if isunix
                [status, ~] = system('which nvidia-smi');
                available = (status == 0);
            else
                available = false;
            end
        end
    end
    
    methods (Access = protected)
        function delete(obj)
            % Destructor - cleanup resources
            if ~isempty(obj.gpuMonitor) && isvalid(obj.gpuMonitor)
                stop(obj.gpuMonitor);
                delete(obj.gpuMonitor);
            end
            
            % Cancel running process
            if obj.isRunning()
                obj.cancel();
            end
        end
    end
end