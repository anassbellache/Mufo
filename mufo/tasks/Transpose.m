classdef Transpose < mufo.core.UFOTask
    % Transpose - UFO transpose task wrapper
    %
    % This task transposes projection data for sinogram generation.
    % Essential for switching between projection and sinogram space.
    %
    % Example:
    %   task = mufo.tasks.Transpose();
    %   task.setDimensions([0, 1, 2]);  % Specify dimension order
    
    properties (Access = private)
        transposeInfo   % Information about transpose operation
    end
    
    methods
        function obj = Transpose(varargin)
            % Constructor
            obj@mufo.core.UFOTask('transpose');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setDimensions(obj, dims)
            % Set dimension ordering for transpose
            %
            % Input:
            %   dims - Array of dimension indices [0, 1, 2] in desired order
            
            if length(dims) ~= 3
                error('Transpose:InvalidDimensions', ...
                    'Dimensions must be a 3-element array');
            end
            
            if ~all(ismember(dims, [0, 1, 2]))
                error('Transpose:InvalidDimensions', ...
                    'Dimensions must contain values 0, 1, 2');
            end
            
            % UFO expects comma-separated string
            dimStr = sprintf('%d,%d,%d', dims(1), dims(2), dims(3));
            obj.setParameter('dimensions', dimStr);
        end
        
        function setFromProjectionToSinogram(obj)
            % Preset for projection to sinogram conversion
            % Swaps projection index with detector row
            obj.setDimensions([1, 0, 2]);
        end
        
        function setFromSinogramToProjection(obj)
            % Preset for sinogram to projection conversion
            % Swaps back to original ordering
            obj.setDimensions([1, 0, 2]);
        end
        
        function info = getTransposeInfo(obj, inputDims)
            % Get information about transpose operation
            %
            % Input:
            %   inputDims - Input data dimensions [x, y, z]
            % Output:
            %   info - Structure with transpose information
            
            dims = obj.getParameter('dimensions');
            if isempty(dims)
                dims = '0,1,2';  % Default no transpose
            end
            
            % Parse dimensions
            dimArray = sscanf(dims, '%d,%d,%d');
            
            info = struct();
            info.inputDims = inputDims;
            info.outputDims = inputDims(dimArray + 1);  % Convert 0-based to 1-based
            info.dimensions = dimArray;
            
            % Determine operation type
            if isequal(dimArray, [0, 1, 2])
                info.operation = 'none';
            elseif isequal(dimArray, [1, 0, 2])
                info.operation = 'projection_to_sinogram';
            elseif isequal(dimArray, [0, 2, 1])
                info.operation = 'swap_y_z';
            else
                info.operation = 'custom';
            end
            
            obj.transposeInfo = info;
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== Transpose Operation ===\n');
                fprintf('Input dimensions: [%d, %d, %d]\n', inputDims);
                fprintf('Output dimensions: [%d, %d, %d]\n', info.outputDims);
                fprintf('Operation: %s\n', info.operation);
                fprintf('Dimension order: [%d, %d, %d]\n', dimArray);
            end
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Check if we need transpose for sinogram generation
                algorithm = config.get('reconstruction.algorithm', 'fbp');
                
                % FBP and similar algorithms need sinogram space
                if any(strcmpi(algorithm, {'fbp', 'sirt', 'sart'}))
                    obj.setFromProjectionToSinogram();
                end
            else
                % Direct struct configuration
                if isfield(config, 'dimensions')
                    if ischar(config.dimensions) || isstring(config.dimensions)
                        obj.setParameter('dimensions', config.dimensions);
                    else
                        obj.setDimensions(config.dimensions);
                    end
                end
                
                if isfield(config, 'operation')
                    switch lower(config.operation)
                        case 'projection_to_sinogram'
                            obj.setFromProjectionToSinogram();
                        case 'sinogram_to_projection'
                            obj.setFromSinogramToProjection();
                    end
                end
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define transpose parameters
            
            % Dimension ordering
            obj.addParameterDef('dimensions', 'string', ...
                'description', 'Dimension ordering (comma-separated)', ...
                'default', '0,1,2');
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            dims = obj.getParameter('dimensions');
            if ~isempty(dims)
                try
                    dimArray = sscanf(dims, '%d,%d,%d');
                    if length(dimArray) ~= 3
                        valid = false;
                        warning('Transpose:InvalidFormat', ...
                            'Dimensions must be three comma-separated integers');
                    elseif ~all(ismember(dimArray, [0, 1, 2]))
                        valid = false;
                        warning('Transpose:InvalidValues', ...
                            'Dimension values must be 0, 1, or 2');
                    elseif length(unique(dimArray)) ~= 3
                        valid = false;
                        warning('Transpose:DuplicateDimensions', ...
                            'Each dimension must appear exactly once');
                    end
                catch
                    valid = false;
                    warning('Transpose:ParseError', ...
                        'Cannot parse dimensions parameter');
                end
            end
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate transpose operations
            
            fprintf('\n=== Transpose Task Demo ===\n\n');
            
            % Create task
            task = mufo.tasks.Transpose();
            
            % Show different transpose operations
            inputDims = [2048, 100, 1800];  % [detector_width, num_projections, detector_height]
            
            fprintf('Original data dimensions: [%d, %d, %d]\n', inputDims);
            fprintf('(detector_width, num_projections, detector_height)\n\n');
            
            % Projection to sinogram
            fprintf('1. Projection to Sinogram transpose:\n');
            task.setFromProjectionToSinogram();
            task.getTransposeInfo(inputDims);
            
            % Custom transpose
            fprintf('\n2. Custom transpose (swap Y and Z):\n');
            task.setDimensions([0, 2, 1]);
            task.getTransposeInfo(inputDims);
            
            % No transpose
            fprintf('\n3. No transpose (identity):\n');
            task.setDimensions([0, 1, 2]);
            task.getTransposeInfo(inputDims);
            
            fprintf('\nDemo complete.\n');
        end
    end
end