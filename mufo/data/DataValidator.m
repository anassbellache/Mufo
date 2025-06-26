classdef DataValidator < handle
    % DataValidator - Data validation and quality assessment for MUFO
    %
    % This class provides validation and quality checks for tomographic data,
    % including consistency checks, outlier detection, and data integrity.
    %
    % Example:
    %   validator = mufo.data.DataValidator();
    %   report = validator.validateProjections(data);
    %   validator.checkFlatFieldQuality(flats, darks);
    
    properties (Access = private)
        lastReport         % Last validation report
        tolerances         % Validation tolerances
        verbosity          % Output verbosity level
    end
    
    properties (Constant)
        DEFAULT_TOLERANCES = struct(...
            'nanFraction', 0.01, ...      % Max 1% NaN values
            'infFraction', 0.001, ...      % Max 0.1% Inf values
            'negativesFraction', 0.01, ... % Max 1% negative values
            'zerosFraction', 0.5, ...      % Max 50% zeros
            'dynamicRange', [10, 1e6], ... % Expected dynamic range
            'snrThreshold', 10, ...        % Minimum SNR
            'flatFieldVariation', 0.2);    % Max flat field variation
    end
    
    events
        ValidationComplete
        ValidationFailed
        WarningIssued
    end
    
    methods
        function obj = DataValidator(varargin)
            % Constructor
            %
            % Parameters:
            %   'Tolerances' - Structure with tolerance values
            %   'Verbosity' - Output level (0-2)
            
            p = inputParser;
            addParameter(p, 'Tolerances', obj.DEFAULT_TOLERANCES, @isstruct);
            addParameter(p, 'Verbosity', 1, @(x) x >= 0 && x <= 2);
            parse(p, varargin{:});
            
            obj.tolerances = obj.mergeTolerances(p.Results.Tolerances);
            obj.verbosity = p.Results.Verbosity;
        end
        
        function report = validateProjections(obj, projections, varargin)
            % Validate projection data
            %
            % Input:
            %   projections - Projection data (3D array)
            %   'Angles' - Projection angles (optional)
            %
            % Output:
            %   report - Validation report structure
            
            p = inputParser;
            addRequired(p, 'projections', @(x) isnumeric(x) && ndims(x) == 3);
            addParameter(p, 'Angles', [], @isnumeric);
            parse(p, projections, varargin{:});
            
            report = struct();
            report.timestamp = datetime('now');
            report.dataSize = size(projections);
            report.dataType = class(projections);
            report.passed = true;
            report.errors = {};
            report.warnings = {};
            report.info = {};
            
            % Basic data checks
            obj.checkDataIntegrity(projections, report);
            
            % Check for common issues
            obj.checkSpecialValues(projections, report);
            obj.checkDynamicRange(projections, report);
            obj.checkNoiseLevel(projections, report);
            
            % Projection-specific checks
            obj.checkProjectionConsistency(projections, report);
            obj.checkForArtifacts(projections, report);
            
            % Angle checks if provided
            if ~isempty(p.Results.Angles)
                obj.checkAngles(p.Results.Angles, size(projections, 1), report);
            end
            
            % Summary
            obj.generateSummary(report);
            
            % Store report
            obj.lastReport = report;
            
            % Notify
            if report.passed
                notify(obj, 'ValidationComplete');
            else
                notify(obj, 'ValidationFailed');
            end
            
            % Display report if verbose
            if obj.verbosity > 0
                obj.displayReport(report);
            end
        end
        
        function report = checkFlatFieldQuality(obj, flats, darks)
            % Check flat field reference quality
            %
            % Inputs:
            %   flats - Flat field images
            %   darks - Dark field images
            %
            % Output:
            %   report - Quality report
            
            report = struct();
            report.timestamp = datetime('now');
            report.passed = true;
            report.errors = {};
            report.warnings = {};
            
            % Check dimensions
            if size(flats, 1) ~= size(darks, 1) || size(flats, 2) ~= size(darks, 2)
                report.errors{end+1} = 'Flat and dark field dimensions do not match';
                report.passed = false;
            end
            
            % Calculate averages
            if ndims(flats) == 3
                avgFlat = mean(flats, 3);
                avgDark = mean(darks, 3);
            else
                avgFlat = flats;
                avgDark = darks;
            end
            
            % Check flat field uniformity
            flatStd = std(avgFlat(:));
            flatMean = mean(avgFlat(:));
            flatCV = flatStd / flatMean;  % Coefficient of variation
            
            report.flatFieldCV = flatCV;
            
            if flatCV > obj.tolerances.flatFieldVariation
                report.warnings{end+1} = sprintf(...
                    'High flat field variation: CV = %.2f (threshold: %.2f)', ...
                    flatCV, obj.tolerances.flatFieldVariation);
            end
            
            % Check dark field stability
            if ndims(darks) == 3
                darkStability = std(mean(mean(darks, 1), 2));
                report.darkStability = darkStability;
                
                if darkStability > 0.1 * mean(avgDark(:))
                    report.warnings{end+1} = 'Dark field shows temporal instability';
                end
            end
            
            % Check signal levels
            signalRange = mean(avgFlat(:)) - mean(avgDark(:));
            report.signalRange = signalRange;
            
            if signalRange < 100
                report.errors{end+1} = 'Insufficient signal range in flat field';
                report.passed = false;
            end
            
            % Check for saturation
            maxVal = obj.getMaxValue(class(flats));
            if any(avgFlat(:) > 0.95 * maxVal)
                report.warnings{end+1} = 'Flat field near saturation';
            end
            
            % Display if verbose
            if obj.verbosity > 0
                obj.displayFlatFieldReport(report);
            end
        end
        
        function report = validateSinogram(obj, sinogram)
            % Validate sinogram data
            %
            % Input:
            %   sinogram - Sinogram data (2D or 3D)
            %
            % Output:
            %   report - Validation report
            
            report = struct();
            report.timestamp = datetime('now');
            report.dataSize = size(sinogram);
            report.passed = true;
            report.errors = {};
            report.warnings = {};
            
            % Basic checks
            obj.checkDataIntegrity(sinogram, report);
            obj.checkSpecialValues(sinogram, report);
            
            % Sinogram-specific checks
            obj.checkSinogramContinuity(sinogram, report);
            obj.checkForRingArtifacts(sinogram, report);
            
            % Display if verbose
            if obj.verbosity > 0
                obj.displayReport(report);
            end
        end
        
        function report = validateReconstruction(obj, volume)
            % Validate reconstructed volume
            %
            % Input:
            %   volume - Reconstructed volume (3D)
            %
            % Output:
            %   report - Validation report
            
            report = struct();
            report.timestamp = datetime('now');
            report.dataSize = size(volume);
            report.passed = true;
            report.errors = {};
            report.warnings = {};
            
            % Basic checks
            obj.checkDataIntegrity(volume, report);
            
            % Reconstruction-specific checks
            obj.checkReconstructionQuality(volume, report);
            obj.checkForReconstructionArtifacts(volume, report);
            
            % Display if verbose
            if obj.verbosity > 0
                obj.displayReport(report);
            end
        end
        
        function compareDatasets(obj, data1, data2, name1, name2)
            % Compare two datasets
            %
            % Inputs:
            %   data1, data2 - Datasets to compare
            %   name1, name2 - Names for display
            
            if nargin < 4
                name1 = 'Dataset 1';
                name2 = 'Dataset 2';
            end
            
            fprintf('\n=== Dataset Comparison ===\n');
            fprintf('Comparing: %s vs %s\n\n', name1, name2);
            
            % Size comparison
            fprintf('Size:\n');
            fprintf('  %s: %s\n', name1, mat2str(size(data1)));
            fprintf('  %s: %s\n', name2, mat2str(size(data2)));
            
            if ~isequal(size(data1), size(data2))
                fprintf('  WARNING: Sizes do not match!\n');
                return;
            end
            
            % Statistics comparison
            stats1 = obj.calculateStatistics(data1);
            stats2 = obj.calculateStatistics(data2);
            
            fprintf('\nStatistics:\n');
            fprintf('  Mean:   %s: %.4f | %s: %.4f | Diff: %.2e\n', ...
                name1, stats1.mean, name2, stats2.mean, ...
                abs(stats1.mean - stats2.mean));
            
            fprintf('  Std:    %s: %.4f | %s: %.4f | Diff: %.2e\n', ...
                name1, stats1.std, name2, stats2.std, ...
                abs(stats1.std - stats2.std));
            
            fprintf('  Min:    %s: %.4f | %s: %.4f\n', ...
                name1, stats1.min, name2, stats2.min);
            
            fprintf('  Max:    %s: %.4f | %s: %.4f\n', ...
                name1, stats1.max, name2, stats2.max);
            
            % Difference metrics
            diff = data1 - data2;
            rmse = sqrt(mean(diff(:).^2));
            nrmse = rmse / (stats1.max - stats1.min);
            maxDiff = max(abs(diff(:)));
            
            fprintf('\nDifference Metrics:\n');
            fprintf('  RMSE: %.4e\n', rmse);
            fprintf('  NRMSE: %.4f%%\n', nrmse * 100);
            fprintf('  Max absolute difference: %.4e\n', maxDiff);
            fprintf('  Correlation: %.6f\n', corr(data1(:), data2(:)));
            
            % Visual comparison if 2D or 3D
            if ndims(data1) >= 2 && obj.verbosity > 1
                obj.visualComparison(data1, data2, name1, name2);
            end
        end
        
        function report = getLastReport(obj)
            % Get last validation report
            report = obj.lastReport;
        end
        
        function saveReport(obj, report, filename)
            % Save validation report
            %
            % Inputs:
            %   report - Validation report
            %   filename - Output filename (.txt, .json, or .html)
            
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case '.txt'
                    obj.saveTextReport(report, filename);
                case '.json'
                    obj.saveJSONReport(report, filename);
                case '.html'
                    obj.saveHTMLReport(report, filename);
                otherwise
                    error('DataValidator:UnsupportedFormat', ...
                        'Unsupported report format: %s', ext);
            end
        end
    end
    
    methods (Access = private)
        function checkDataIntegrity(obj, data, report)
            % Basic data integrity checks
            
            % Check if empty
            if isempty(data)
                report.errors{end+1} = 'Data is empty';
                report.passed = false;
                return;
            end
            
            % Check dimensions
            if any(size(data) == 0)
                report.errors{end+1} = 'Data has zero dimension';
                report.passed = false;
            end
            
            % Check data type
            if ~isnumeric(data)
                report.errors{end+1} = 'Data is not numeric';
                report.passed = false;
            end
        end
        
        function checkSpecialValues(obj, data, report)
            % Check for special values (NaN, Inf, negatives)
            
            totalElements = numel(data);
            
            % NaN values
            nanCount = sum(isnan(data(:)));
            nanFraction = nanCount / totalElements;
            report.nanCount = nanCount;
            report.nanFraction = nanFraction;
            
            if nanFraction > obj.tolerances.nanFraction
                report.errors{end+1} = sprintf(...
                    'Too many NaN values: %.2f%% (threshold: %.2f%%)', ...
                    nanFraction * 100, obj.tolerances.nanFraction * 100);
                report.passed = false;
            elseif nanCount > 0
                report.warnings{end+1} = sprintf(...
                    'Data contains %d NaN values (%.4f%%)', ...
                    nanCount, nanFraction * 100);
            end
            
            % Inf values
            infCount = sum(isinf(data(:)));
            infFraction = infCount / totalElements;
            report.infCount = infCount;
            report.infFraction = infFraction;
            
            if infFraction > obj.tolerances.infFraction
                report.errors{end+1} = sprintf(...
                    'Too many Inf values: %.2f%% (threshold: %.2f%%)', ...
                    infFraction * 100, obj.tolerances.infFraction * 100);
                report.passed = false;
            elseif infCount > 0
                report.warnings{end+1} = sprintf(...
                    'Data contains %d Inf values (%.4f%%)', ...
                    infCount, infFraction * 100);
            end
            
            % Negative values (for projection data)
            negCount = sum(data(:) < 0);
            negFraction = negCount / totalElements;
            report.negativeCount = negCount;
            report.negativeFraction = negFraction;
            
            if negFraction > obj.tolerances.negativesFraction
                report.warnings{end+1} = sprintf(...
                    'Many negative values: %.2f%% (threshold: %.2f%%)', ...
                    negFraction * 100, obj.tolerances.negativesFraction * 100);
            end
            
            % Zero values
            zeroCount = sum(data(:) == 0);
            zeroFraction = zeroCount / totalElements;
            report.zeroCount = zeroCount;
            report.zeroFraction = zeroFraction;
            
            if zeroFraction > obj.tolerances.zerosFraction
                report.warnings{end+1} = sprintf(...
                    'Many zero values: %.2f%% (threshold: %.2f%%)', ...
                    zeroFraction * 100, obj.tolerances.zerosFraction * 100);
            end
        end
        
        function checkDynamicRange(obj, data, report)
            % Check data dynamic range
            
            validData = data(~isnan(data(:)) & ~isinf(data(:)));
            
            if isempty(validData)
                report.errors{end+1} = 'No valid data for dynamic range check';
                report.passed = false;
                return;
            end
            
            minVal = min(validData);
            maxVal = max(validData);
            dynamicRange = maxVal / max(minVal, eps);
            
            report.dynamicRange = dynamicRange;
            report.dataRange = [minVal, maxVal];
            
            if dynamicRange < obj.tolerances.dynamicRange(1)
                report.warnings{end+1} = sprintf(...
                    'Low dynamic range: %.1f (expected > %.1f)', ...
                    dynamicRange, obj.tolerances.dynamicRange(1));
            elseif dynamicRange > obj.tolerances.dynamicRange(2)
                report.warnings{end+1} = sprintf(...
                    'Very high dynamic range: %.1e (expected < %.1e)', ...
                    dynamicRange, obj.tolerances.dynamicRange(2));
            end
        end
        
        function checkNoiseLevel(obj, data, report)
            % Estimate and check noise level
            
            % Estimate noise from high-frequency components
            if ndims(data) == 3
                % Use difference between adjacent projections
                if size(data, 1) > 1
                    diffs = diff(data, 1, 1);
                    noise = std(diffs(:)) / sqrt(2);
                else
                    % Use single slice
                    slice = data(:, :, round(size(data, 3)/2));
                    noise = obj.estimateNoise2D(slice);
                end
            else
                noise = obj.estimateNoise2D(data);
            end
            
            % Calculate SNR
            signal = mean(abs(data(:)));
            snr = signal / noise;
            
            report.noiseLevel = noise;
            report.snr = snr;
            
            if snr < obj.tolerances.snrThreshold
                report.warnings{end+1} = sprintf(...
                    'Low SNR: %.1f (threshold: %.1f)', ...
                    snr, obj.tolerances.snrThreshold);
            end
        end
        
        function checkProjectionConsistency(obj, projections, report)
            % Check consistency across projections
            
            [ny, nx, nz] = size(projections);
            
            % Check mean intensity variations
            meanIntensities = squeeze(mean(mean(projections, 1), 2));
            intensityCV = std(meanIntensities) / mean(meanIntensities);
            
            report.intensityCV = intensityCV;
            
            if intensityCV > 0.1
                report.warnings{end+1} = sprintf(...
                    'High intensity variation across projections: CV = %.2f', ...
                    intensityCV);
            end
            
            % Check for sudden jumps
            diffs = abs(diff(meanIntensities));
            threshold = 3 * median(diffs);
            jumps = find(diffs > threshold);
            
            if ~isempty(jumps)
                report.warnings{end+1} = sprintf(...
                    'Detected %d intensity jumps at projections: %s', ...
                    length(jumps), mat2str(jumps));
            end
        end
        
        function checkForArtifacts(obj, data, report)
            % Check for common artifacts
            
            % Check for dead pixels (constant across projections)
            if ndims(data) == 3
                pixelStd = std(data, 0, 1);
                deadPixels = sum(pixelStd(:) == 0);
                
                if deadPixels > 0
                    report.warnings{end+1} = sprintf(...
                        'Detected %d dead pixels', deadPixels);
                    report.deadPixelMap = (pixelStd == 0);
                end
                
                % Check for hot pixels
                medianStd = median(pixelStd(:));
                hotPixels = sum(pixelStd(:) > 10 * medianStd);
                
                if hotPixels > 0
                    report.warnings{end+1} = sprintf(...
                        'Detected %d hot pixels', hotPixels);
                end
            end
        end
        
        function checkAngles(obj, angles, nProjections, report)
            % Check projection angles
            
            if length(angles) ~= nProjections
                report.errors{end+1} = sprintf(...
                    'Number of angles (%d) does not match projections (%d)', ...
                    length(angles), nProjections);
                report.passed = false;
                return;
            end
            
            % Check angle range
            angleRange = max(angles) - min(angles);
            report.angleRange = angleRange;
            
            % Check spacing
            angleSteps = diff(angles);
            
            if std(angleSteps) > 0.01 * mean(angleSteps)
                report.warnings{end+1} = 'Non-uniform angle spacing detected';
            end
            
            % Check for duplicate angles
            uniqueAngles = length(unique(angles));
            if uniqueAngles < length(angles)
                report.warnings{end+1} = sprintf(...
                    'Duplicate angles detected: %d unique out of %d', ...
                    uniqueAngles, length(angles));
            end
        end
        
        function checkSinogramContinuity(obj, sinogram, report)
            % Check sinogram continuity
            
            if size(sinogram, 1) < 2
                return;
            end
            
            % Check first and last projections (should be similar for 360°)
            if ndims(sinogram) == 2
                firstProj = sinogram(1, :);
                lastProj = sinogram(end, :);
            else
                firstProj = squeeze(sinogram(1, :, round(size(sinogram, 3)/2)));
                lastProj = squeeze(sinogram(end, :, round(size(sinogram, 3)/2)));
            end
            
            correlation = corr(firstProj(:), lastProj(:));
            report.endCorrelation = correlation;
            
            if correlation < 0.9
                report.warnings{end+1} = sprintf(...
                    'Low correlation between first/last projection: %.3f', ...
                    correlation);
            end
        end
        
        function checkForRingArtifacts(obj, sinogram, report)
            % Check for ring artifacts in sinogram
            
            % Analyze vertical lines in sinogram
            if ndims(sinogram) == 2
                verticalMean = mean(sinogram, 1);
                verticalStd = std(sinogram, 0, 1);
            else
                % Use middle slice
                slice = sinogram(:, :, round(size(sinogram, 3)/2));
                verticalMean = mean(slice, 1);
                verticalStd = std(slice, 0, 1);
            end
            
            % Look for outliers in std (rings have low std)
            medianStd = median(verticalStd);
            ringCandidates = find(verticalStd < 0.1 * medianStd);
            
            if length(ringCandidates) > 5
                report.warnings{end+1} = sprintf(...
                    'Potential ring artifacts detected at %d positions', ...
                    length(ringCandidates));
                report.ringPositions = ringCandidates;
            end
        end
        
        function checkReconstructionQuality(obj, volume, report)
            % Check reconstruction quality metrics
            
            % Check value range
            stats = obj.calculateStatistics(volume);
            report.stats = stats;
            
            % Check for reconstruction artifacts
            % Look for negative values in regions that should be positive
            negativeVolume = sum(volume(:) < 0) / numel(volume);
            
            if negativeVolume > 0.1
                report.warnings{end+1} = sprintf(...
                    'Large negative region in reconstruction: %.1f%%', ...
                    negativeVolume * 100);
            end
            
            % Check sharpness (using gradient magnitude)
            midSlice = volume(:, :, round(size(volume, 3)/2));
            [Gmag, ~] = imgradient(midSlice);
            sharpness = mean(Gmag(:));
            report.sharpness = sharpness;
        end
        
        function checkForReconstructionArtifacts(obj, volume, report)
            % Check for common reconstruction artifacts
            
            % Check for cupping artifacts (higher values at edges)
            midSlice = volume(:, :, round(size(volume, 3)/2));
            center = midSlice(round(end/2), round(end/2));
            edges = [midSlice(1, :), midSlice(end, :), ...
                    midSlice(:, 1)', midSlice(:, end)'];
            edgeMean = mean(edges);
            
            if edgeMean > 1.2 * center
                report.warnings{end+1} = 'Possible cupping artifact detected';
            end
            
            % Check for streak artifacts
            % Simplified check using variance along radial lines
            [ny, nx] = size(midSlice);
            cy = ny / 2;
            cx = nx / 2;
            
            radialVars = zeros(8, 1);
            angles = linspace(0, pi, 8);
            
            for i = 1:length(angles)
                % Extract radial profile
                [x, y] = pol2cart(angles(i), 0:min(cx, cy));
                x = round(x + cx);
                y = round(y + cy);
                
                valid = x >= 1 & x <= nx & y >= 1 & y <= ny;
                profile = midSlice(sub2ind([ny, nx], y(valid), x(valid)));
                
                radialVars(i) = var(profile);
            end
            
            if max(radialVars) / min(radialVars) > 10
                report.warnings{end+1} = 'Possible streak artifacts detected';
            end
        end
        
        function noise = estimateNoise2D(obj, image)
            % Estimate noise level in 2D image
            
            % Use median absolute deviation of high-pass filtered image
            h = [-1 2 -1; 2 -4 2; -1 2 -1] / 4;  % Laplacian
            filtered = conv2(double(image), h, 'valid');
            
            % Robust noise estimate
            noise = 1.4826 * median(abs(filtered(:) - median(filtered(:))));
        end
        
        function stats = calculateStatistics(obj, data)
            % Calculate comprehensive statistics
            
            validData = data(~isnan(data(:)) & ~isinf(data(:)));
            
            stats = struct();
            stats.mean = mean(validData);
            stats.std = std(validData);
            stats.min = min(validData);
            stats.max = max(validData);
            stats.median = median(validData);
            stats.p01 = prctile(validData, 1);
            stats.p99 = prctile(validData, 99);
            stats.iqr = iqr(validData);
            stats.skewness = skewness(validData);
            stats.kurtosis = kurtosis(validData);
        end
        
        function maxVal = getMaxValue(obj, dataType)
            % Get maximum value for data type
            
            switch dataType
                case 'uint8'
                    maxVal = 255;
                case 'uint16'
                    maxVal = 65535;
                case 'uint32'
                    maxVal = 4294967295;
                case 'int8'
                    maxVal = 127;
                case 'int16'
                    maxVal = 32767;
                case 'int32'
                    maxVal = 2147483647;
                otherwise
                    maxVal = inf;
            end
        end
        
        function merged = mergeTolerances(obj, custom)
            % Merge custom tolerances with defaults
            
            merged = obj.DEFAULT_TOLERANCES;
            fields = fieldnames(custom);
            
            for i = 1:length(fields)
                if isfield(merged, fields{i})
                    merged.(fields{i}) = custom.(fields{i});
                end
            end
        end
        
        function generateSummary(obj, report)
            % Generate summary section of report
            
            report.summary = struct();
            report.summary.numErrors = length(report.errors);
            report.summary.numWarnings = length(report.warnings);
            
            if report.passed
                if report.summary.numWarnings == 0
                    report.summary.status = 'PASSED - No issues found';
                else
                    report.summary.status = sprintf(...
                        'PASSED with %d warnings', report.summary.numWarnings);
                end
            else
                report.summary.status = sprintf(...
                    'FAILED - %d errors, %d warnings', ...
                    report.summary.numErrors, report.summary.numWarnings);
            end
        end
        
        function displayReport(obj, report)
            % Display validation report
            
            fprintf('\n=== Data Validation Report ===\n');
            fprintf('Timestamp: %s\n', datestr(report.timestamp));
            fprintf('Data size: %s\n', mat2str(report.dataSize));
            fprintf('Data type: %s\n', report.dataType);
            fprintf('\nStatus: %s\n', report.summary.status);
            
            if ~isempty(report.errors)
                fprintf('\nERRORS:\n');
                for i = 1:length(report.errors)
                    fprintf('  - %s\n', report.errors{i});
                end
            end
            
            if ~isempty(report.warnings)
                fprintf('\nWARNINGS:\n');
                for i = 1:length(report.warnings)
                    fprintf('  - %s\n', report.warnings{i});
                end
            end
            
            if obj.verbosity > 1
                fprintf('\nDETAILS:\n');
                
                if isfield(report, 'nanCount')
                    fprintf('  NaN values: %d (%.4f%%)\n', ...
                        report.nanCount, report.nanFraction * 100);
                end
                
                if isfield(report, 'dynamicRange')
                    fprintf('  Dynamic range: %.1e\n', report.dynamicRange);
                    fprintf('  Data range: [%.3e, %.3e]\n', ...
                        report.dataRange(1), report.dataRange(2));
                end
                
                if isfield(report, 'snr')
                    fprintf('  SNR: %.1f\n', report.snr);
                    fprintf('  Noise level: %.3e\n', report.noiseLevel);
                end
            end
            
            fprintf('\n');
        end
        
        function displayFlatFieldReport(obj, report)
            % Display flat field quality report
            
            fprintf('\n=== Flat Field Quality Report ===\n');
            fprintf('Timestamp: %s\n', datestr(report.timestamp));
            fprintf('Status: %s\n', report.passed);
            
            fprintf('\nMETRICS:\n');
            fprintf('  Flat field CV: %.3f\n', report.flatFieldCV);
            fprintf('  Signal range: %.1f\n', report.signalRange);
            
            if isfield(report, 'darkStability')
                fprintf('  Dark stability: %.3f\n', report.darkStability);
            end
            
            if ~isempty(report.errors)
                fprintf('\nERRORS:\n');
                for i = 1:length(report.errors)
                    fprintf('  - %s\n', report.errors{i});
                end
            end
            
            if ~isempty(report.warnings)
                fprintf('\nWARNINGS:\n');
                for i = 1:length(report.warnings)
                    fprintf('  - %s\n', report.warnings{i});
                end
            end
            
            fprintf('\n');
        end
        
        function visualComparison(obj, data1, data2, name1, name2)
            % Visual comparison of datasets
            
            figure('Name', 'Dataset Comparison');
            
            if ndims(data1) == 3
                % Use middle slice
                slice1 = data1(:, :, round(size(data1, 3)/2));
                slice2 = data2(:, :, round(size(data2, 3)/2));
            else
                slice1 = data1;
                slice2 = data2;
            end
            
            % Show images
            subplot(2,3,1);
            imagesc(slice1);
            axis image;
            colorbar;
            title(name1);
            
            subplot(2,3,2);
            imagesc(slice2);
            axis image;
            colorbar;
            title(name2);
            
            subplot(2,3,3);
            imagesc(slice1 - slice2);
            axis image;
            colorbar;
            title('Difference');
            
            % Histograms
            subplot(2,3,4);
            histogram(slice1(:), 100);
            title([name1, ' Histogram']);
            
            subplot(2,3,5);
            histogram(slice2(:), 100);
            title([name2, ' Histogram']);
            
            subplot(2,3,6);
            histogram((slice1(:) - slice2(:)), 100);
            title('Difference Histogram');
            
            colormap gray;
        end
        
        function saveTextReport(obj, report, filename)
            % Save report as text file
            
            fid = fopen(filename, 'w');
            if fid == -1
                error('DataValidator:FileError', ...
                    'Cannot create file: %s', filename);
            end
            
            fprintf(fid, 'Data Validation Report\n');
            fprintf(fid, '======================\n\n');
            fprintf(fid, 'Generated: %s\n', datestr(report.timestamp));
            fprintf(fid, 'Data size: %s\n', mat2str(report.dataSize));
            fprintf(fid, 'Data type: %s\n\n', report.dataType);
            
            fprintf(fid, 'STATUS: %s\n\n', report.summary.status);
            
            if ~isempty(report.errors)
                fprintf(fid, 'ERRORS:\n');
                for i = 1:length(report.errors)
                    fprintf(fid, '  * %s\n', report.errors{i});
                end
                fprintf(fid, '\n');
            end
            
            if ~isempty(report.warnings)
                fprintf(fid, 'WARNINGS:\n');
                for i = 1:length(report.warnings)
                    fprintf(fid, '  * %s\n', report.warnings{i});
                end
                fprintf(fid, '\n');
            end
            
            % Write all numeric fields
            fprintf(fid, 'METRICS:\n');
            fields = fieldnames(report);
            for i = 1:length(fields)
                if isnumeric(report.(fields{i})) && isscalar(report.(fields{i}))
                    fprintf(fid, '  %s: %g\n', fields{i}, report.(fields{i}));
                end
            end
            
            fclose(fid);
        end
        
        function saveJSONReport(obj, report, filename)
            % Save report as JSON file
            
            % Clean report for JSON (remove non-serializable fields)
            cleanReport = report;
            fields = fieldnames(cleanReport);
            
            for i = 1:length(fields)
                if ~(isnumeric(cleanReport.(fields{i})) || ...
                     ischar(cleanReport.(fields{i})) || ...
                     islogical(cleanReport.(fields{i})) || ...
                     iscell(cleanReport.(fields{i})) || ...
                     isstruct(cleanReport.(fields{i})))
                    cleanReport = rmfield(cleanReport, fields{i});
                end
            end
            
            % Convert datetime
            if isfield(cleanReport, 'timestamp')
                cleanReport.timestamp = datestr(cleanReport.timestamp);
            end
            
            % Write JSON
            jsonStr = jsonencode(cleanReport, 'PrettyPrint', true);
            
            fid = fopen(filename, 'w');
            fprintf(fid, '%s', jsonStr);
            fclose(fid);
        end
        
        function saveHTMLReport(obj, report, filename)
            % Save report as HTML file
            
            fid = fopen(filename, 'w');
            if fid == -1
                error('DataValidator:FileError', ...
                    'Cannot create file: %s', filename);
            end
            
            % HTML header
            fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n');
            fprintf(fid, '<title>Data Validation Report</title>\n');
            fprintf(fid, '<style>\n');
            fprintf(fid, 'body { font-family: Arial, sans-serif; margin: 20px; }\n');
            fprintf(fid, 'h1, h2 { color: #333; }\n');
            fprintf(fid, '.error { color: #d00; }\n');
            fprintf(fid, '.warning { color: #f80; }\n');
            fprintf(fid, '.passed { color: #080; }\n');
            fprintf(fid, '.failed { color: #d00; }\n');
            fprintf(fid, 'table { border-collapse: collapse; margin: 10px 0; }\n');
            fprintf(fid, 'th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n');
            fprintf(fid, 'th { background-color: #f5f5f5; }\n');
            fprintf(fid, '</style>\n</head>\n<body>\n');
            
            % Content
            fprintf(fid, '<h1>Data Validation Report</h1>\n');
            fprintf(fid, '<p><strong>Generated:</strong> %s</p>\n', datestr(report.timestamp));
            fprintf(fid, '<p><strong>Data size:</strong> %s</p>\n', mat2str(report.dataSize));
            fprintf(fid, '<p><strong>Data type:</strong> %s</p>\n', report.dataType);
            
            % Status
            if report.passed
                statusClass = 'passed';
            else
                statusClass = 'failed';
            end
            fprintf(fid, '<h2>Status: <span class="%s">%s</span></h2>\n', ...
                statusClass, report.summary.status);
            
            % Errors
            if ~isempty(report.errors)
                fprintf(fid, '<h2>Errors</h2>\n<ul>\n');
                for i = 1:length(report.errors)
                    fprintf(fid, '<li class="error">%s</li>\n', report.errors{i});
                end
                fprintf(fid, '</ul>\n');
            end
            
            % Warnings
            if ~isempty(report.warnings)
                fprintf(fid, '<h2>Warnings</h2>\n<ul>\n');
                for i = 1:length(report.warnings)
                    fprintf(fid, '<li class="warning">%s</li>\n', report.warnings{i});
                end
                fprintf(fid, '</ul>\n');
            end
            
            % Metrics table
            fprintf(fid, '<h2>Metrics</h2>\n');
            fprintf(fid, '<table>\n<tr><th>Metric</th><th>Value</th></tr>\n');
            
            fields = fieldnames(report);
            for i = 1:length(fields)
                if isnumeric(report.(fields{i})) && isscalar(report.(fields{i}))
                    fprintf(fid, '<tr><td>%s</td><td>%.4g</td></tr>\n', ...
                        fields{i}, report.(fields{i}));
                end
            end
            
            fprintf(fid, '</table>\n');
            
            % HTML footer
            fprintf(fid, '</body>\n</html>\n');
            fclose(fid);
        end
    end
    
    methods (Static)
        function demo()
            % Demonstrate DataValidator capabilities
            
            fprintf('\n=== DataValidator Demo ===\n\n');
            
            % Create validator
            validator = mufo.data.DataValidator('Verbosity', 2);
            
            % Create test data with various issues
            testData = phantom(256) * 1000;
            testData = repmat(testData, [1, 1, 50]);
            
            % Add some issues
            testData(50:60, 50:60, :) = NaN;  % NaN region
            testData(100, 100, :) = 0;        % Dead pixel
            testData(150, 150, :) = 5000;     % Hot pixel
            testData = testData + 10 * randn(size(testData));  % Noise
            
            % Example 1: Validate projections
            fprintf('1. Validating projection data:\n');
            report = validator.validateProjections(testData);
            
            % Example 2: Check flat field quality
            fprintf('\n2. Checking flat field quality:\n');
            flats = ones(256, 256, 10) * 1000 + 50 * randn(256, 256, 10);
            darks = ones(256, 256, 10) * 100 + 10 * randn(256, 256, 10);
            
            flatReport = validator.checkFlatFieldQuality(flats, darks);
            
            % Example 3: Compare datasets
            fprintf('\n3. Comparing datasets:\n');
            testData2 = testData + 5 * randn(size(testData));
            validator.compareDatasets(testData, testData2, 'Original', 'Noisy');
            
            fprintf('\nDemo complete.\n');
        end
    end
end