classdef Filter < mufo.core.UFOTask
    % Filter - UFO filter task wrapper
    %
    % Applies various filters for tomographic reconstruction.
    % Typically used in Fourier space between FFT and IFFT operations.
    %
    % Example:
    %   task = mufo.tasks.Filter();
    %   task.setFilter('butterworth');
    %   task.setCutoff(0.5);
    %   task.setOrder(4);
    
    properties (Access = private)
        filterInfo      % Information about filter configuration
    end
    
    methods
        function obj = Filter(varargin)
            % Constructor
            obj@mufo.core.UFOTask('filter');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setFilter(obj, filterType)
            % Set filter type
            %
            % Input:
            %   filterType - Filter name ('ramp', 'butterworth', 'hamming', etc.)
            
            validFilters = {'ramp', 'ramp-fromreal', 'butterworth', ...
                           'hamming', 'hann', 'cosine', 'shepp-logan'};
            
            if ~any(strcmpi(filterType, validFilters))
                warning('Filter:UnknownType', ...
                    'Unknown filter type: %s', filterType);
            end
            
            obj.setParameter('filter', filterType);
        end
        
        function setCutoff(obj, cutoff)
            % Set cutoff frequency (normalized 0-1)
            %
            % Input:
            %   cutoff - Cutoff frequency (0-1, where 1 is Nyquist)
            
            if cutoff <= 0 || cutoff > 1
                error('Filter:InvalidCutoff', ...
                    'Cutoff must be between 0 and 1');
            end
            
            obj.setParameter('cutoff', cutoff);
        end
        
        function setOrder(obj, order)
            % Set filter order (for Butterworth)
            %
            % Input:
            %   order - Filter order (positive integer)
            
            if order < 1 || mod(order, 1) ~= 0
                error('Filter:InvalidOrder', ...
                    'Order must be a positive integer');
            end
            
            obj.setParameter('order', order);
        end
        
        function setSigma(obj, sigma)
            % Set sigma parameter (for Gaussian-based filters)
            %
            % Input:
            %   sigma - Gaussian sigma value
            
            if sigma <= 0
                error('Filter:InvalidSigma', ...
                    'Sigma must be positive');
            end
            
            obj.setParameter('sigma', sigma);
        end
        
        function [H, f] = generateFilterResponse(obj, N)
            % Generate filter frequency response
            %
            % Input:
            %   N - Number of frequency points
            % Outputs:
            %   H - Filter response
            %   f - Normalized frequency axis (0-1)
            
            if nargin < 2
                N = 512;
            end
            
            filterType = obj.getParameter('filter');
            if isempty(filterType)
                filterType = 'ramp-fromreal';
            end
            
            cutoff = obj.getParameter('cutoff');
            if isempty(cutoff)
                cutoff = 1.0;
            end
            
            % Create frequency axis (0 to 1)
            f = linspace(0, 1, N/2);
            
            % Generate filter response
            switch lower(filterType)
                case 'ramp'
                    H = f;
                    
                case 'ramp-fromreal'
                    H = f;
                    H(1) = 0;  % Zero DC component
                    
                case 'butterworth'
                    order = obj.getParameter('order');
                    if isempty(order), order = 4; end
                    H = 1 ./ (1 + (f/cutoff).^(2*order));
                    H = H .* f;  % Combine with ramp
                    
                case 'hamming'
                    H = f .* (0.54 + 0.46 * cos(pi * f / cutoff));
                    H(f > cutoff) = 0;
                    
                case 'hann'
                    H = f .* 0.5 * (1 + cos(pi * f / cutoff));
                    H(f > cutoff) = 0;
                    
                case 'cosine'
                    H = f .* cos(pi * f / (2 * cutoff));
                    H(f > cutoff) = 0;
                    
                case 'shepp-logan'
                    % Shepp-Logan filter: |f| * sinc(f/(2*fc))
                    H = f .* abs(sinc(f / (2 * cutoff)));
                    
                otherwise
                    H = f;  % Default to ramp
                    warning('Filter:UnknownType', ...
                        'Unknown filter type, using ramp');
            end
            
            % Store filter info
            obj.filterInfo = struct(...
                'type', filterType, ...
                'cutoff', cutoff, ...
                'response', H, ...
                'frequency', f);
        end
        
        function previewFilter(obj, N)
            % Preview filter frequency response
            %
            % Input:
            %   N - Number of frequency points (default: 1024)
            
            if nargin < 2
                N = 1024;
            end
            
            % Generate filter response
            [H, f] = obj.generateFilterResponse(N);
            
            % Create figure
            figure('Name', 'Filter Preview');
            
            % Plot frequency response
            subplot(2,1,1);
            plot(f, H, 'b-', 'LineWidth', 2);
            xlabel('Normalized Frequency');
            ylabel('Filter Response');
            title(sprintf('%s Filter', obj.getParameter('filter')));
            grid on;
            
            % Add cutoff line
            cutoff = obj.getParameter('cutoff');
            if ~isempty(cutoff) && cutoff < 1
                hold on;
                plot([cutoff, cutoff], [0, max(H)], 'r--', 'LineWidth', 1);
                legend('Filter', sprintf('Cutoff = %.2f', cutoff), ...
                    'Location', 'best');
            end
            
            % Plot in dB
            subplot(2,1,2);
            HdB = 20*log10(abs(H + eps));  % Add eps to avoid log(0)
            plot(f, HdB, 'b-', 'LineWidth', 2);
            xlabel('Normalized Frequency');
            ylabel('Filter Response (dB)');
            grid on;
            ylim([-60, max(HdB)+5]);
            
            % Add cutoff line
            if ~isempty(cutoff) && cutoff < 1
                hold on;
                plot([cutoff, cutoff], ylim, 'r--', 'LineWidth', 1);
            end
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get reconstruction filter settings
                filterType = config.get('reconstruction.filter', 'ramp-fromreal');
                obj.setFilter(filterType);
                
                % Get filter parameters
                cutoff = config.get('reconstruction.filterCutoff', 1.0);
                if cutoff < 1
                    obj.setCutoff(cutoff);
                end
                
                % Butterworth order
                if strcmpi(filterType, 'butterworth')
                    order = config.get('reconstruction.filterOrder', 4);
                    obj.setOrder(order);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'filter')
                    obj.setFilter(config.filter);
                end
                
                if isfield(config, 'cutoff')
                    obj.setCutoff(config.cutoff);
                end
                
                if isfield(config, 'order')
                    obj.setOrder(config.order);
                end
                
                if isfield(config, 'sigma')
                    obj.setSigma(config.sigma);
                end
            end
        end
        
        function compareFilters(obj, filters, N)
            % Compare multiple filter responses
            %
            % Inputs:
            %   filters - Cell array of filter configurations
            %   N - Number of frequency points
            
            if nargin < 3
                N = 1024;
            end
            
            figure('Name', 'Filter Comparison');
            colors = lines(length(filters));
            
            % Store original settings
            originalParams = obj.parameters;
            
            hold on;
            legendEntries = {};
            
            for i = 1:length(filters)
                % Apply filter configuration
                filterConfig = filters{i};
                
                if isfield(filterConfig, 'type')
                    obj.setFilter(filterConfig.type);
                end
                if isfield(filterConfig, 'cutoff')
                    obj.setCutoff(filterConfig.cutoff);
                end
                if isfield(filterConfig, 'order')
                    obj.setOrder(filterConfig.order);
                end
                
                % Generate response
                [H, f] = obj.generateFilterResponse(N);
                
                % Plot
                plot(f, H, 'Color', colors(i,:), 'LineWidth', 2);
                
                % Legend entry
                legendEntries{i} = sprintf('%s (cutoff=%.2f)', ...
                    obj.getParameter('filter'), ...
                    obj.getParameter('cutoff'));
            end
            
            % Restore original settings
            obj.parameters = originalParams;
            
            % Finalize plot
            xlabel('Normalized Frequency');
            ylabel('Filter Response');
            title('Filter Comparison');
            legend(legendEntries, 'Location', 'best');
            grid on;
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define filter parameters
            
            % Filter type
            obj.addParameterDef('filter', 'string', ...
                'description', 'Filter type', ...
                'values', {'ramp', 'ramp-fromreal', 'butterworth', ...
                          'hamming', 'hann', 'cosine', 'shepp-logan'}, ...
                'default', 'ramp-fromreal');
            
            % Cutoff frequency
            obj.addParameterDef('cutoff', 'numeric', ...
                'description', 'Cutoff frequency (0-1)', ...
                'range', [0, 1], ...
                'default', 1.0);
            
            % Filter order (for Butterworth)
            obj.addParameterDef('order', 'integer', ...
                'description', 'Filter order (Butterworth)', ...
                'range', [1, 10], ...
                'default', 4);
            
            % Sigma (for Gaussian-based filters)
            obj.addParameterDef('sigma', 'positive', ...
                'description', 'Gaussian sigma parameter', ...
                'default', 1.0);
            
            % Padding mode
            obj.addParameterDef('padding-mode', 'string', ...
                'description', 'Padding mode for filtering', ...
                'values', {'none', 'zero', 'replicate'}, ...
                'default', 'zero');
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            filterType = obj.getParameter('filter');
            if isempty(filterType)
                % Filter type will default to ramp-fromreal
                return;
            end
            
            % Check filter-specific parameters
            switch lower(filterType)
                case 'butterworth'
                    order = obj.getParameter('order');
                    if isempty(order)
                        warning('Filter:MissingOrder', ...
                            'Butterworth filter requires order parameter');
                    end
                    
                case {'hamming', 'hann', 'cosine', 'shepp-logan'}
                    cutoff = obj.getParameter('cutoff');
                    if isempty(cutoff) || cutoff >= 1
                        warning('Filter:NoCutoff', ...
                            '%s filter typically uses cutoff < 1', filterType);
                    end
            end
        end
    end
    
    methods (Access = private)
        function y = sinc(obj, x)
            % Sinc function implementation
            y = ones(size(x));
            idx = x ~= 0;
            y(idx) = sin(pi * x(idx)) ./ (pi * x(idx));
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate filter operations
            
            fprintf('\n=== Filter Task Demo ===\n\n');
            
            % Create task
            task = mufo.tasks.Filter();
            
            % Example 1: Standard FBP filters
            fprintf('1. Previewing standard FBP filter (ramp-fromreal):\n');
            task.setFilter('ramp-fromreal');
            task.previewFilter();
            
            % Example 2: Butterworth filter
            fprintf('\n2. Butterworth filter with different orders:\n');
            filters = {
                struct('type', 'butterworth', 'cutoff', 0.5, 'order', 2), ...
                struct('type', 'butterworth', 'cutoff', 0.5, 'order', 4), ...
                struct('type', 'butterworth', 'cutoff', 0.5, 'order', 8)
            };
            task.compareFilters(filters);
            
            % Example 3: Different filter types
            fprintf('\n3. Comparing different filter types:\n');
            filters = {
                struct('type', 'ramp-fromreal', 'cutoff', 1.0), ...
                struct('type', 'hamming', 'cutoff', 0.7), ...
                struct('type', 'hann', 'cutoff', 0.7), ...
                struct('type', 'shepp-logan', 'cutoff', 0.8)
            };
            task.compareFilters(filters);
            
            fprintf('\nDemo complete.\n');
        end
    end
end