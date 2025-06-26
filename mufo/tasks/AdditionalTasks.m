% This file contains multiple task classes for common UFO operations

%% Transpose Task
classdef Transpose < mufo.core.UFOTask
    % Transpose - UFO transpose task wrapper
    %
    % Transposes projection data for sinogram generation
    
    methods
        function obj = Transpose(varargin)
            obj@mufo.core.UFOTask('transpose');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % No specific parameters for transpose
        end
    end
end

%% FFT Task
classdef FFT < mufo.core.UFOTask
    % FFT - UFO FFT task wrapper
    %
    % Performs Fast Fourier Transform along specified dimensions
    
    methods
        function obj = FFT(varargin)
            obj@mufo.core.UFOTask('fft');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setDimensions(obj, dims)
            % Set FFT dimensions (1, 2, or 3)
            obj.setParameter('dimensions', dims);
        end
        
        function setSize(obj, sizes)
            % Set FFT sizes for each dimension
            if length(sizes) == 2
                obj.setParameter('size-x', sizes(1));
                obj.setParameter('size-y', sizes(2));
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            obj.addParameterDef('dimensions', 'integer', ...
                'description', 'Number of FFT dimensions', ...
                'values', [1, 2, 3], ...
                'default', 1);
            
            obj.addParameterDef('size-x', 'integer', ...
                'description', 'FFT size in X', ...
                'range', [1, inf]);
            
            obj.addParameterDef('size-y', 'integer', ...
                'description', 'FFT size in Y', ...
                'range', [1, inf]);
        end
    end
end

%% IFFT Task
classdef IFFT < mufo.core.UFOTask
    % IFFT - UFO inverse FFT task wrapper
    
    methods
        function obj = IFFT(varargin)
            obj@mufo.core.UFOTask('ifft');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setDimensions(obj, dims)
            obj.setParameter('dimensions', dims);
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            obj.addParameterDef('dimensions', 'integer', ...
                'description', 'Number of IFFT dimensions', ...
                'values', [1, 2, 3], ...
                'default', 1);
        end
    end
end

%% Filter Task
classdef Filter < mufo.core.UFOTask
    % Filter - UFO filter task wrapper
    %
    % Applies various filters for tomographic reconstruction
    
    methods
        function obj = Filter(varargin)
            obj@mufo.core.UFOTask('filter');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setFilter(obj, filterType)
            % Set filter type
            obj.setParameter('filter', filterType);
        end
        
        function setCutoff(obj, cutoff)
            % Set cutoff frequency (0-1)
            obj.setParameter('cutoff', cutoff);
        end
        
        function setOrder(obj, order)
            % Set filter order (for Butterworth)
            obj.setParameter('order', order);
        end
        
        function previewFilter(obj, size)
            % Preview filter frequency response
            if nargin < 2
                size = 512;
            end
            
            filterType = obj.getParameter('filter');
            cutoff = obj.getParameter('cutoff');
            
            % Create frequency axis
            f = linspace(0, 1, size/2);
            
            % Generate filter response
            switch filterType
                case 'ramp'
                    H = f;
                case 'ramp-fromreal'
                    H = f;
                    H(1) = 0;
                case 'butterworth'
                    order = obj.getParameter('order');
                    if isempty(order), order = 4; end
                    if isempty(cutoff), cutoff = 0.5; end
                    H = 1 ./ (1 + (f/cutoff).^(2*order));
                    H = H .* f;  % Combine with ramp
                case 'hamming'
                    if isempty(cutoff), cutoff = 0.5; end
                    H = f .* (0.54 + 0.46 * cos(pi * f / cutoff));
                    H(f > cutoff) = 0;
                case 'hann'
                    if isempty(cutoff), cutoff = 0.5; end
                    H = f .* 0.5 * (1 + cos(pi * f / cutoff));
                    H(f > cutoff) = 0;
                otherwise
                    H = f;  % Default ramp
            end
            
            % Plot
            figure('Name', 'Filter Preview');
            plot(f, H, 'b-', 'LineWidth', 2);
            xlabel('Normalized Frequency');
            ylabel('Filter Response');
            title(sprintf('%s Filter', filterType));
            grid on;
            
            if ~isempty(cutoff)
                hold on;
                plot([cutoff, cutoff], [0, max(H)], 'r--', 'LineWidth', 1);
                legend('Filter', sprintf('Cutoff = %.2f', cutoff));
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            obj.addParameterDef('filter', 'string', ...
                'description', 'Filter type', ...
                'values', {'ramp', 'ramp-fromreal', 'butterworth', ...
                          'hamming', 'hann', 'cosine', 'shepp-logan'}, ...
                'default', 'ramp-fromreal');
            
            obj.addParameterDef('cutoff', 'numeric', ...
                'description', 'Cutoff frequency (0-1)', ...
                'range', [0, 1], ...
                'default', 1.0);
            
            obj.addParameterDef('order', 'integer', ...
                'description', 'Filter order (Butterworth)', ...
                'range', [1, 10], ...
                'default', 4);
        end
    end
end

%% Backproject Task
classdef Backproject < mufo.core.UFOTask
    % Backproject - UFO backprojection task wrapper
    %
    % Performs filtered backprojection reconstruction
    
    methods
        function obj = Backproject(varargin)
            obj@mufo.core.UFOTask('backproject');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setCenter(obj, center)
            % Set center of rotation
            obj.setParameter('axis-pos', center);
        end
        
        function setAngleStep(obj, step)
            % Set angle step in radians
            obj.setParameter('angle-step', step);
        end
        
        function setAngleRange(obj, startAngle, endAngle)
            % Set angle range
            obj.setParameter('angle-start', startAngle);
            obj.setParameter('angle-end', endAngle);
        end
        
        function setPrecision(obj, mode)
            % Set precision mode
            obj.setParameter('precision-mode', mode);
        end
        
        function configureFromJSON(obj, config)
            % Configure from JSON config
            if isa(config, 'mufo.Config')
                reconParams = config.getReconstructionParams();
                
                if ~isempty(reconParams.center)
                    obj.setCenter(reconParams.center);
                end
                
                % Calculate angle step
                angleRange = reconParams.angleRange;
                if isfield(reconParams, 'numProjections')
                    angleStep = deg2rad(angleRange / reconParams.numProjections);
                    obj.setAngleStep(angleStep);
                end
                
                if isfield(reconParams, 'precision')
                    obj.setPrecision(reconParams.precision);
                end
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            obj.addParameterDef('axis-pos', 'numeric', ...
                'description', 'Axis position (center of rotation)', ...
                'required', true);
            
            obj.addParameterDef('angle-step', 'numeric', ...
                'description', 'Angle step in radians', ...
                'default', []);
            
            obj.addParameterDef('angle-start', 'numeric', ...
                'description', 'Starting angle in radians', ...
                'default', 0);
            
            obj.addParameterDef('angle-end', 'numeric', ...
                'description', 'Ending angle in radians', ...
                'default', 2*pi);
            
            obj.addParameterDef('precision-mode', 'string', ...
                'description', 'Precision mode', ...
                'values', {'single', 'half', 'int8'}, ...
                'default', 'single');
            
            obj.addParameterDef('addressing-mode', 'string', ...
                'description', 'Addressing mode', ...
                'values', {'clamp', 'clamp_to_edge', 'repeat', 'mirrored_repeat'}, ...
                'default', 'clamp_to_edge');
        end
        
        function valid = validate(obj)
            valid = true;
            
            % Check if center is set
            if isempty(obj.getParameter('axis-pos'))
                valid = false;
                warning('Backproject:NoCenter', ...
                    'Center of rotation must be set');
            end
            
            % Check angle step
            angleStep = obj.getParameter('angle-step');
            if isempty(angleStep)
                % Try to calculate from angle range
                angleStart = obj.getParameter('angle-start');
                angleEnd = obj.getParameter('angle-end');
                
                warning('Backproject:NoAngleStep', ...
                    'Angle step not set - UFO will use default');
            end
        end
    end
end

%% Ring Remove Task
classdef RingRemove < mufo.core.UFOTask
    % RingRemove - UFO ring removal task wrapper
    
    methods
        function obj = RingRemove(varargin)
            obj@mufo.core.UFOTask('remove-rings');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setMethod(obj, method)
            % Set ring removal method
            obj.setParameter('method', method);
        end
        
        function setThreshold(obj, threshold)
            % Set threshold for ring detection
            obj.setParameter('threshold', threshold);
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            obj.addParameterDef('method', 'string', ...
                'description', 'Ring removal method', ...
                'values', {'wavelet', 'morphological'}, ...
                'default', 'wavelet');
            
            obj.addParameterDef('threshold', 'numeric', ...
                'description', 'Detection threshold', ...
                'range', [0, inf], ...
                'default', 1.0);
            
            obj.addParameterDef('sigma', 'numeric', ...
                'description', 'Gaussian sigma', ...
                'range', [0, inf], ...
                'default', 1.0);
        end
    end
end

%% Average Task
classdef Average < mufo.core.UFOTask
    % Average - UFO averaging task wrapper
    %
    % Averages multiple frames (e.g., for dark/flat fields)
    
    methods
        function obj = Average(varargin)
            obj@mufo.core.UFOTask('average');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % No specific parameters - averages all input
        end
    end
end

%% Stack Task
classdef Stack < mufo.core.UFOTask
    % Stack - UFO stack task wrapper
    %
    % Stacks multiple inputs into single output
    
    methods
        function obj = Stack(varargin)
            obj@mufo.core.UFOTask('stack');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setNumber(obj, num)
            % Set number of items to stack
            obj.setParameter('number', num);
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            obj.addParameterDef('number', 'integer', ...
                'description', 'Number of items to stack', ...
                'range', [1, inf], ...
                'default', 1);
        end
    end
end

%% Crop Task
classdef Crop < mufo.core.UFOTask
    % Crop - UFO crop task wrapper
    
    methods
        function obj = Crop(varargin)
            obj@mufo.core.UFOTask('crop');
            obj.defineParameters();
            
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setRegion(obj, x, y, width, height)
            % Set crop region
            obj.setParameter('x', x);
            obj.setParameter('y', y);
            obj.setParameter('width', width);
            obj.setParameter('height', height);
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            obj.addParameterDef('x', 'integer', ...
                'description', 'X offset', ...
                'range', [0, inf], ...
                'required', true);
            
            obj.addParameterDef('y', 'integer', ...
                'description', 'Y offset', ...
                'range', [0, inf], ...
                'required', true);
            
            obj.addParameterDef('width', 'integer', ...
                'description', 'Crop width', ...
                'range', [1, inf], ...
                'required', true);
            
            obj.addParameterDef('height', 'integer', ...
                'description', 'Crop height', ...
                'range', [1, inf], ...
                'required', true);
        end
    end
end