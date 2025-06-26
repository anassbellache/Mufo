classdef PhaseRetrieve < mufo.core.UFOTask
    % PhaseRetrieve - UFO Paganin phase retrieval task wrapper
    %
    % This task performs phase retrieval using the Paganin method,
    % which is particularly effective for inline phase contrast imaging.
    %
    % Example:
    %   task = mufo.tasks.PhaseRetrieve();
    %   task.setEnergy(25);          % 25 keV
    %   task.setDistance(1.5);       % 1.5 meters
    %   task.setPixelSize(1.1e-6);   % 1.1 microns
    
    properties (Access = private)
        materialDB      % Material properties database
        wavelength      % Calculated wavelength
    end
    
    properties (Constant)
        % Physical constants
        PLANCK_EV = 4.135667696e-15;  % Planck constant in eV·s
        SPEED_OF_LIGHT = 299792458;    % Speed of light in m/s
    end
    
    methods
        function obj = PhaseRetrieve(varargin)
            % Constructor
            obj@mufo.core.UFOTask('retrieve-phase');
            
            % Define parameters
            obj.defineParameters();
            
            % Initialize material database
            obj.initializeMaterialDB();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setEnergy(obj, energy)
            % Set X-ray energy in keV
            %
            % Input:
            %   energy - X-ray energy in keV
            
            if energy <= 0
                error('PhaseRetrieve:InvalidEnergy', ...
                    'Energy must be positive');
            end
            
            obj.setParameter('energy', energy);
            
            % Calculate wavelength
            obj.wavelength = obj.calculateWavelength(energy);
        end
        
        function setDistance(obj, distance)
            % Set propagation distance in meters
            %
            % Input:
            %   distance - Sample-to-detector distance in meters
            
            if distance <= 0
                error('PhaseRetrieve:InvalidDistance', ...
                    'Distance must be positive');
            end
            
            obj.setParameter('distance', distance);
        end
        
        function setPixelSize(obj, pixelSize)
            % Set detector pixel size in meters
            %
            % Input:
            %   pixelSize - Pixel size in meters
            
            if pixelSize <= 0
                error('PhaseRetrieve:InvalidPixelSize', ...
                    'Pixel size must be positive');
            end
            
            obj.setParameter('pixel-size', pixelSize);
        end
        
        function setMaterial(obj, material)
            % Set material for automatic delta/beta calculation
            %
            % Input:
            %   material - Material name or formula (e.g., 'H2O', 'Al')
            
            if ischar(material) || isstring(material)
                if isfield(obj.materialDB, upper(material))
                    matData = obj.materialDB.(upper(material));
                    
                    % Calculate delta and beta for current energy
                    energy = obj.getParameter('energy');
                    if ~isempty(energy)
                        [delta, beta] = obj.calculateOpticalConstants(...
                            material, energy);
                        
                        obj.setParameter('delta', delta);
                        obj.setParameter('beta', beta);
                        
                        fprintf('Material %s at %.1f keV:\n', ...
                            material, energy);
                        fprintf('  delta = %.3e\n', delta);
                        fprintf('  beta = %.3e\n', beta);
                    else
                        warning('PhaseRetrieve:NoEnergy', ...
                            'Set energy before material for automatic calculation');
                    end
                else
                    warning('PhaseRetrieve:UnknownMaterial', ...
                        'Unknown material: %s', material);
                end
            end
        end
        
        function setOpticalConstants(obj, delta, beta)
            % Manually set optical constants
            %
            % Inputs:
            %   delta - Real part decrement of refractive index
            %   beta - Imaginary part of refractive index
            
            obj.setParameter('delta', delta);
            obj.setParameter('beta', beta);
        end
        
        function setRegularization(obj, method, parameter)
            % Set regularization method and parameter
            %
            % Inputs:
            %   method - 'TIE' or 'CTF'
            %   parameter - Regularization parameter
            
            obj.setParameter('regularization-mode', method);
            obj.setParameter('regularization-parameter', parameter);
        end
        
        function setThickness(obj, thickness)
            % Set sample thickness for quantitative retrieval
            %
            % Input:
            %   thickness - Sample thickness in meters
            
            obj.setParameter('thickness', thickness);
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get Paganin parameters
                params = config.getPaganinParams();
                
                if isfield(params, 'enabled') && ~params.enabled
                    warning('PhaseRetrieve:Disabled', ...
                        'Phase retrieval is disabled in configuration');
                end
                
                obj.setEnergy(params.energy);
                obj.setDistance(params.distance);
                obj.setPixelSize(params.pixelSize);
                
                % Set optical constants
                if isfield(params, 'material') && ~isempty(params.material)
                    obj.setMaterial(params.material);
                else
                    obj.setOpticalConstants(params.delta, params.beta);
                end
                
                % Apply log setting
                obj.setParameter('apply-log', params.applyLog);
                
            else
                % Direct struct configuration
                if isfield(config, 'energy')
                    obj.setEnergy(config.energy);
                end
                
                if isfield(config, 'distance')
                    obj.setDistance(config.distance);
                end
                
                if isfield(config, 'pixelSize')
                    obj.setPixelSize(config.pixelSize);
                end
                
                if isfield(config, 'delta')
                    obj.setParameter('delta', config.delta);
                end
                
                if isfield(config, 'beta')
                    obj.setParameter('beta', config.beta);
                end
                
                if isfield(config, 'applyLog')
                    obj.setParameter('apply-log', config.applyLog);
                end
            end
        end
        
        function params = calculatePaganinParameters(obj)
            % Calculate Paganin filter parameters
            %
            % Output:
            %   params - Structure with calculated parameters
            
            energy = obj.getParameter('energy');
            distance = obj.getParameter('distance');
            pixelSize = obj.getParameter('pixel-size');
            delta = obj.getParameter('delta');
            beta = obj.getParameter('beta');
            
            % Check all required parameters
            if isempty(energy) || isempty(distance) || isempty(pixelSize)
                error('PhaseRetrieve:MissingParameters', ...
                    'Energy, distance, and pixel size are required');
            end
            
            if isempty(delta) || isempty(beta)
                error('PhaseRetrieve:MissingOpticalConstants', ...
                    'Delta and beta are required');
            end
            
            % Calculate wavelength if not done
            if isempty(obj.wavelength)
                obj.wavelength = obj.calculateWavelength(energy);
            end
            
            % Calculate Paganin parameter
            % α = (δ/β) * λ * z / (4π)
            alpha = (delta / beta) * obj.wavelength * distance / (4 * pi);
            
            % Calculate cutoff frequency
            cutoffFreq = sqrt(beta / (delta * obj.wavelength * distance)) / pixelSize;
            
            % Create output structure
            params = struct();
            params.energy = energy;
            params.wavelength = obj.wavelength;
            params.distance = distance;
            params.pixelSize = pixelSize;
            params.delta = delta;
            params.beta = beta;
            params.alpha = alpha;
            params.cutoffFrequency = cutoffFreq;
            params.fresnel = obj.wavelength * distance / (pixelSize^2);
            
            % Display if no output requested
            if nargout == 0
                fprintf('\n=== Paganin Parameters ===\n');
                fprintf('Energy: %.1f keV\n', energy);
                fprintf('Wavelength: %.3e m\n', obj.wavelength);
                fprintf('Distance: %.3f m\n', distance);
                fprintf('Pixel size: %.2f µm\n', pixelSize * 1e6);
                fprintf('Delta: %.3e\n', delta);
                fprintf('Beta: %.3e\n', beta);
                fprintf('Alpha parameter: %.3e\n', alpha);
                fprintf('Cutoff frequency: %.3f /pixel\n', cutoffFreq);
                fprintf('Fresnel number: %.3f\n', params.fresnel);
            end
        end
        
        function previewPhaseRetrieval(obj, inputData)
            % Preview phase retrieval on sample data
            %
            % Input:
            %   inputData - 2D input image
            
            % Get parameters
            params = obj.calculatePaganinParameters();
            
            % Apply Paganin filter (simplified MATLAB version)
            [ny, nx] = size(inputData);
            
            % Create frequency grid
            fx = ifftshift((-nx/2:nx/2-1) / nx);
            fy = ifftshift((-ny/2:ny/2-1) / ny);
            [FX, FY] = meshgrid(fx, fy);
            
            % Frequency squared
            f2 = FX.^2 + FY.^2;
            
            % Paganin filter
            % H = 1 / (1 + α * λ * z * f²)
            H = 1 ./ (1 + params.alpha * params.wavelength * params.distance * ...
                     f2 / (params.pixelSize^2));
            
            % Apply filter in Fourier domain
            fftData = fft2(inputData);
            filtered = real(ifft2(fftData .* H));
            
            % Apply logarithm if requested
            if obj.getParameter('apply-log')
                filtered = -log(filtered);
                filtered(filtered < 0) = 0;
            end
            
            % Display comparison
            figure('Name', 'Phase Retrieval Preview');
            
            subplot(2,3,1);
            imagesc(inputData);
            colormap gray;
            axis image;
            colorbar;
            title('Input');
            
            subplot(2,3,2);
            imagesc(filtered);
            colormap gray;
            axis image;
            colorbar;
            title('Phase Retrieved');
            
            subplot(2,3,3);
            imagesc(filtered - inputData);
            colormap gray;
            axis image;
            colorbar;
            title('Difference');
            
            % Show frequency filter
            subplot(2,3,4);
            imagesc(fftshift(H));
            colormap hot;
            axis image;
            colorbar;
            title('Paganin Filter');
            
            % Show profiles
            subplot(2,3,5);
            midLine = round(ny/2);
            plot(inputData(midLine, :), 'b-', 'DisplayName', 'Input');
            hold on;
            plot(filtered(midLine, :), 'r-', 'DisplayName', 'Retrieved');
            xlabel('Pixel');
            ylabel('Intensity');
            legend();
            title('Horizontal Profile');
            
            % Show parameter info
            subplot(2,3,6);
            axis off;
            text(0.1, 0.9, sprintf('Energy: %.1f keV', params.energy));
            text(0.1, 0.8, sprintf('Distance: %.2f m', params.distance));
            text(0.1, 0.7, sprintf('Pixel: %.1f µm', params.pixelSize*1e6));
            text(0.1, 0.6, sprintf('δ/β: %.1f', params.delta/params.beta));
            text(0.1, 0.5, sprintf('α: %.2e', params.alpha));
            title('Parameters');
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define phase retrieval parameters
            
            % X-ray parameters
            obj.addParameterDef('energy', 'positive', ...
                'description', 'X-ray energy in keV', ...
                'required', true);
            
            % Geometry parameters
            obj.addParameterDef('distance', 'positive', ...
                'description', 'Propagation distance in meters', ...
                'required', true);
            
            obj.addParameterDef('pixel-size', 'positive', ...
                'description', 'Pixel size in meters', ...
                'required', true);
            
            % Material properties
            obj.addParameterDef('delta', 'positive', ...
                'description', 'Real part decrement of refractive index', ...
                'required', true);
            
            obj.addParameterDef('beta', 'positive', ...
                'description', 'Imaginary part of refractive index', ...
                'required', true);
            
            % Optional parameters
            obj.addParameterDef('apply-log', 'logical', ...
                'description', 'Apply logarithm after retrieval', ...
                'default', true);
            
            obj.addParameterDef('regularization-mode', 'string', ...
                'description', 'Regularization method', ...
                'values', {'TIE', 'CTF'}, ...
                'default', 'TIE');
            
            obj.addParameterDef('regularization-parameter', 'positive', ...
                'description', 'Regularization parameter', ...
                'default', 0.5);
            
            obj.addParameterDef('thickness', 'positive', ...
                'description', 'Sample thickness in meters');
            
            % Padding
            obj.addParameterDef('padding', 'integer', ...
                'description', 'Padding mode', ...
                'values', [0, 1, 2], ...
                'default', 1);
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            % Check required parameters
            energy = obj.getParameter('energy');
            distance = obj.getParameter('distance');
            pixelSize = obj.getParameter('pixel-size');
            delta = obj.getParameter('delta');
            beta = obj.getParameter('beta');
            
            if isempty(energy) || isempty(distance) || isempty(pixelSize)
                valid = false;
                return;
            end
            
            if isempty(delta) || isempty(beta)
                valid = false;
                warning('PhaseRetrieve:MissingOpticalConstants', ...
                    'Delta and beta must be set');
                return;
            end
            
            % Check delta/beta ratio
            if delta/beta < 10
                warning('PhaseRetrieve:LowContrast', ...
                    'Low delta/beta ratio (%.1f) may result in poor phase retrieval', ...
                    delta/beta);
            end
            
            % Check Fresnel number
            wavelength = obj.calculateWavelength(energy);
            fresnel = wavelength * distance / (pixelSize^2);
            
            if fresnel < 0.1
                warning('PhaseRetrieve:LowFresnel', ...
                    'Low Fresnel number (%.3f) - phase effects may be weak', ...
                    fresnel);
            elseif fresnel > 10
                warning('PhaseRetrieve:HighFresnel', ...
                    'High Fresnel number (%.3f) - consider near-field approximation', ...
                    fresnel);
            end
        end
    end
    
    methods (Access = private)
        function wavelength = calculateWavelength(obj, energy_keV)
            % Calculate X-ray wavelength from energy
            %
            % Input:
            %   energy_keV - Energy in keV
            % Output:
            %   wavelength - Wavelength in meters
            
            energy_J = energy_keV * 1000 * 1.602176634e-19;  % Convert to Joules
            wavelength = obj.PLANCK_EV * 1000 * obj.SPEED_OF_LIGHT / (energy_keV * 1000);
        end
        
        function initializeMaterialDB(obj)
            % Initialize material properties database
            obj.materialDB = struct();
            
            % Common materials (simplified - real values depend on energy)
            % Format: density (g/cm³), Z_effective
            obj.materialDB.H2O = struct('density', 1.0, 'Z', 7.42);
            obj.materialDB.AL = struct('density', 2.70, 'Z', 13);
            obj.materialDB.SI = struct('density', 2.33, 'Z', 14);
            obj.materialDB.FE = struct('density', 7.87, 'Z', 26);
            obj.materialDB.CU = struct('density', 8.96, 'Z', 29);
            obj.materialDB.BONE = struct('density', 1.92, 'Z', 13.8);
            obj.materialDB.TISSUE = struct('density', 1.06, 'Z', 7.6);
            obj.materialDB.PMMA = struct('density', 1.19, 'Z', 6.56);
        end
        
        function [delta, beta] = calculateOpticalConstants(obj, material, energy_keV)
            % Calculate optical constants for material at given energy
            % This is a simplified calculation - use tabulated values for accuracy
            
            if ~isfield(obj.materialDB, upper(material))
                error('PhaseRetrieve:UnknownMaterial', ...
                    'Unknown material: %s', material);
            end
            
            mat = obj.materialDB.(upper(material));
            
            % Classical electron radius
            r_e = 2.818e-15;  % meters
            
            % Wavelength
            lambda = obj.calculateWavelength(energy_keV);
            
            % Number density (atoms/m³)
            N_A = 6.022e23;  % Avogadro's number
            n = mat.density * 1e6 * N_A / (mat.Z * 1.008);  % Simplified
            
            % Delta (simplified - ignores anomalous dispersion)
            delta = n * r_e * lambda^2 / (2 * pi);
            
            % Beta (simplified - uses empirical scaling)
            % β ∝ λ⁴ for photoelectric absorption
            beta = delta / (1000 * (energy_keV / 10)^3);  % Very rough approximation
            
            % Scale by typical values
            delta = delta * 1e6;  % Typical values are ~1e-6
            beta = beta * 1e3;    % Typical values are ~1e-9
        end
    end
    
    methods (Static)
        function demo()
            % Run a demonstration of phase retrieval
            
            fprintf('\n=== Phase Retrieval Demo ===\n\n');
            
            % Create task
            task = mufo.tasks.PhaseRetrieve();
            
            % Set typical parameters
            task.setEnergy(25);           % 25 keV
            task.setDistance(1.5);        % 1.5 m
            task.setPixelSize(6.5e-6);    % 6.5 µm
            
            % Set material
            task.setMaterial('PMMA');
            
            % Show parameters
            task.calculatePaganinParameters();
            
            % Create synthetic phase contrast image
            fprintf('\nCreating synthetic phase contrast data...\n');
            
            % Create phantom
            [X, Y] = meshgrid(-256:255, -256:255);
            R = sqrt(X.^2 + Y.^2);
            
            % Phase object (circles with different phases)
            phase = zeros(512, 512);
            phase(R < 100) = 0.5;
            phase(R < 50) = 1.0;
            phase(abs(X) < 20 & abs(Y) < 150) = 0.3;
            
            % Add edge enhancement (phase contrast)
            edgeEnhanced = phase - 0.1 * del2(phase);
            edgeEnhanced = edgeEnhanced + 1;  % Add background
            
            % Add noise
            edgeEnhanced = edgeEnhanced + 0.02 * randn(512, 512);
            
            % Preview retrieval
            fprintf('Previewing phase retrieval...\n');
            task.previewPhaseRetrieval(edgeEnhanced);
            
            fprintf('\nDemo complete.\n');
        end
    end
end