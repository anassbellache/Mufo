classdef FlatFieldCorrect < mufo.core.UFOTask
    % FlatFieldCorrect - UFO flat-field correction task wrapper
    %
    % This task performs flat-field correction using dark and flat
    % reference images. It requires exactly 3 inputs: projections,
    % darks, and flats.
    %
    % Example:
    %   task = mufo.tasks.FlatFieldCorrect();
    %   task.setAbsorptionCorrect(true);
    %   task.setDarkScale(0.95);
    
    properties (Access = private)
        correctionInfo  % Information about correction parameters
    end
    
    methods
        function obj = FlatFieldCorrect(varargin)
            % Constructor
            obj@mufo.core.UFOTask('flat-field-correct');
            
            % Define parameters
            obj.defineParameters();
            
            % Set parameters if provided
            if ~isempty(varargin)
                obj.setParameters(varargin{:});
            end
        end
        
        function setAbsorptionCorrect(obj, enable)
            % Enable/disable absorption correction
            %
            % Input:
            %   enable - true for -log((P-D)/(F-D)), false for (P-D)/(F-D)
            
            obj.setParameter('absorption-correct', enable);
        end
        
        function setFixNanInf(obj, enable)
            % Enable/disable NaN and Inf fixing
            %
            % Input:
            %   enable - Fix invalid values (default: true)
            
            obj.setParameter('fix-nan-and-inf', enable);
        end
        
        function setDarkScale(obj, scale)
            % Set dark field scaling factor
            %
            % Input:
            %   scale - Scaling factor for dark field (default: 1.0)
            
            if scale <= 0
                error('FlatFieldCorrect:InvalidScale', ...
                    'Dark scale must be positive');
            end
            
            obj.setParameter('dark-scale', scale);
        end
        
        function setFlatScale(obj, scale)
            % Set flat field scaling factor
            %
            % Input:
            %   scale - Scaling factor for flat field (default: 1.0)
            
            if scale <= 0
                error('FlatFieldCorrect:InvalidScale', ...
                    'Flat scale must be positive');
            end
            
            obj.setParameter('flat-scale', scale);
        end
        
        function setNanValue(obj, value)
            % Set replacement value for NaN
            %
            % Input:
            %   value - Value to replace NaN with
            
            obj.setParameter('nan-value', value);
        end
        
        function setInfValue(obj, value)
            % Set replacement value for Inf
            %
            % Input:
            %   value - Value to replace Inf with
            
            obj.setParameter('inf-value', value);
        end
        
        function configureFromJSON(obj, config)
            % Configure task from JSON config
            %
            % Input:
            %   config - Config object or struct
            
            if isa(config, 'mufo.Config')
                % Get flat-field settings
                ffConfig = config.get('preprocessing.flatField', struct());
                
                if isfield(ffConfig, 'absorptionCorrect')
                    obj.setAbsorptionCorrect(ffConfig.absorptionCorrect);
                end
                
                if isfield(ffConfig, 'fixNanAndInf')
                    obj.setFixNanInf(ffConfig.fixNanAndInf);
                end
                
                if isfield(ffConfig, 'darkScale')
                    obj.setDarkScale(ffConfig.darkScale);
                end
                
                if isfield(ffConfig, 'flatScale')
                    obj.setFlatScale(ffConfig.flatScale);
                end
                
            else
                % Direct struct configuration
                if isfield(config, 'absorptionCorrect')
                    obj.setAbsorptionCorrect(config.absorptionCorrect);
                end
                
                if isfield(config, 'fixNanAndInf')
                    obj.setFixNanInf(config.fixNanAndInf);
                end
                
                if isfield(config, 'darkScale')
                    obj.setDarkScale(config.darkScale);
                end
                
                if isfield(config, 'flatScale')
                    obj.setFlatScale(config.flatScale);
                end
                
                if isfield(config, 'nanValue')
                    obj.setNanValue(config.nanValue);
                end
                
                if isfield(config, 'infValue')
                    obj.setInfValue(config.infValue);
                end
            end
        end
        
        function chain = createCorrectionChain(obj, projPath, darkPath, flatPath, outputPath)
            % Create complete flat-field correction chain
            %
            % Inputs:
            %   projPath - Path to projections
            %   darkPath - Path to dark fields
            %   flatPath - Path to flat fields
            %   outputPath - Output path
            % Output:
            %   chain - UFOChain object
            
            chain = mufo.core.UFOChain('Name', 'FlatFieldCorrection');
            
            % Create read tasks
            projRead = mufo.tasks.Read();
            projRead.setPath(projPath);
            
            darkRead = mufo.tasks.Read();
            darkRead.setPath(darkPath);
            
            flatRead = mufo.tasks.Read();
            flatRead.setPath(flatPath);
            
            % Add parallel read tasks
            chain.addParallel({projRead, darkRead, flatRead});
            
            % Add flat-field correction
            chain.addTask(obj);
            
            % Add write task
            writeTask = mufo.tasks.Write();
            writeTask.setFilename(outputPath);
            chain.addTask(writeTask);
        end
        
        function validateInputs(obj, projInfo, darkInfo, flatInfo)
            % Validate input data compatibility
            %
            % Inputs:
            %   projInfo - Projection data info
            %   darkInfo - Dark field info
            %   flatInfo - Flat field info
            
            % Check dimensions compatibility
            if darkInfo.dims(1) ~= projInfo.dims(1) || ...
               darkInfo.dims(2) ~= projInfo.dims(2)
                warning('FlatFieldCorrect:DimensionMismatch', ...
                    'Dark field dimensions do not match projections');
            end
            
            if flatInfo.dims(1) ~= projInfo.dims(1) || ...
               flatInfo.dims(2) ~= projInfo.dims(2)
                warning('FlatFieldCorrect:DimensionMismatch', ...
                    'Flat field dimensions do not match projections');
            end
            
            % Store correction info
            obj.correctionInfo = struct();
            obj.correctionInfo.projDims = projInfo.dims;
            obj.correctionInfo.numDarks = darkInfo.dims(3);
            obj.correctionInfo.numFlats = flatInfo.dims(3);
            
            % Display info
            fprintf('\nFlat-field correction setup:\n');
            fprintf('  Projections: %dx%dx%d\n', projInfo.dims);
            fprintf('  Dark fields: %d frames\n', obj.correctionInfo.numDarks);
            fprintf('  Flat fields: %d frames\n', obj.correctionInfo.numFlats);
            fprintf('  Absorption correction: %s\n', ...
                mat2str(obj.getParameter('absorption-correct')));
        end
        
        function previewCorrection(obj, projData, darkData, flatData)
            % Preview flat-field correction on sample data
            %
            % Inputs:
            %   projData - Sample projection (2D)
            %   darkData - Dark field (2D)
            %   flatData - Flat field (2D)
            
            % Get parameters
            absorptionCorrect = obj.getParameter('absorption-correct');
            darkScale = obj.getParameter('dark-scale');
            flatScale = obj.getParameter('flat-scale');
            
            % Apply correction
            if absorptionCorrect
                % Absorption correction: -log((P-D)/(F-D))
                corrected = -log((projData - darkScale*darkData) ./ ...
                                (flatScale*flatData - darkScale*darkData));
            else
                % Simple correction: (P-D)/(F-D)
                corrected = (projData - darkScale*darkData) ./ ...
                           (flatScale*flatData - darkScale*darkData);
            end
            
            % Fix NaN and Inf if enabled
            if obj.getParameter('fix-nan-and-inf')
                nanValue = obj.getParameter('nan-value');
                infValue = obj.getParameter('inf-value');
                
                if isempty(nanValue)
                    nanValue = 0;
                end
                if isempty(infValue)
                    infValue = 0;
                end
                
                corrected(isnan(corrected)) = nanValue;
                corrected(isinf(corrected)) = infValue;
            end
            
            % Display comparison
            figure('Name', 'Flat-field Correction Preview');
            
            subplot(2,2,1);
            imagesc(projData);
            colormap gray;
            axis image;
            colorbar;
            title('Original Projection');
            
            subplot(2,2,2);
            imagesc(darkData);
            colormap gray;
            axis image;
            colorbar;
            title('Dark Field');
            
            subplot(2,2,3);
            imagesc(flatData);
            colormap gray;
            axis image;
            colorbar;
            title('Flat Field');
            
            subplot(2,2,4);
            imagesc(corrected);
            colormap gray;
            axis image;
            colorbar;
            title('Corrected');
            
            % Display statistics
            fprintf('\nCorrection statistics:\n');
            fprintf('  Original range: [%.3f, %.3f]\n', ...
                min(projData(:)), max(projData(:)));
            fprintf('  Corrected range: [%.3f, %.3f]\n', ...
                min(corrected(:)), max(corrected(:)));
            fprintf('  NaN pixels: %d\n', sum(isnan(corrected(:))));
            fprintf('  Inf pixels: %d\n', sum(isinf(corrected(:))));
        end
    end
    
    methods (Access = protected)
        function defineParameters(obj)
            % Define flat-field correction parameters
            
            % Absorption correction
            obj.addParameterDef('absorption-correct', 'logical', ...
                'description', 'Apply absorption correction (-log)', ...
                'default', true);
            
            % Fix invalid values
            obj.addParameterDef('fix-nan-and-inf', 'logical', ...
                'description', 'Fix NaN and Inf values', ...
                'default', true);
            
            % Scaling factors
            obj.addParameterDef('dark-scale', 'numeric', ...
                'description', 'Dark field scaling factor', ...
                'default', 1.0, ...
                'range', [0, inf]);
            
            obj.addParameterDef('flat-scale', 'numeric', ...
                'description', 'Flat field scaling factor', ...
                'default', 1.0, ...
                'range', [0, inf]);
            
            % Replacement values
            obj.addParameterDef('nan-value', 'numeric', ...
                'description', 'Value to replace NaN', ...
                'default', 0);
            
            obj.addParameterDef('inf-value', 'numeric', ...
                'description', 'Value to replace Inf', ...
                'default', 0);
            
            % Error handling
            obj.addParameterDef('error-mode', 'string', ...
                'description', 'Error handling mode', ...
                'values', {'ignore', 'warn', 'error'}, ...
                'default', 'warn');
        end
        
        function valid = validate(obj)
            % Custom validation
            valid = true;
            
            % Check scaling factors
            darkScale = obj.getParameter('dark-scale');
            flatScale = obj.getParameter('flat-scale');
            
            if darkScale <= 0 || flatScale <= 0
                valid = false;
                warning('FlatFieldCorrect:InvalidScale', ...
                    'Scaling factors must be positive');
            end
            
            % Note: This task requires exactly 3 inputs (proj, dark, flat)
            % but input validation is handled by UFO at runtime
        end
    end
    
    methods (Static)
        function demo()
            % Run a demonstration of flat-field correction
            
            fprintf('\n=== Flat-field Correction Demo ===\n\n');
            
            % Create synthetic data
            fprintf('Creating synthetic data...\n');
            
            % Projection with some structure
            [X, Y] = meshgrid(1:512, 1:512);
            proj = 1000 + 500*exp(-((X-256).^2 + (Y-256).^2)/10000);
            
            % Add some noise
            proj = proj + 50*randn(512, 512);
            
            % Dark field (offset)
            dark = 100 + 10*randn(512, 512);
            
            % Flat field (gain variations)
            flat = 2000 + 200*sin(X/50) + 200*cos(Y/50) + 50*randn(512, 512);
            
            % Create task and preview
            task = mufo.tasks.FlatFieldCorrect();
            task.setAbsorptionCorrect(true);
            task.setFixNanInf(true);
            
            fprintf('Previewing correction...\n');
            task.previewCorrection(proj, dark, flat);
            
            fprintf('\nDemo complete.\n');
        end
    end
end