% test_hdf5_ffc_pipeline.m
% Test pipeline for HDF5 data loading, flat field correction, and visualization
% Tests the complete workflow from HDF5 input to corrected output with inspection

%% Setup and Configuration
clear; close all;
fprintf('MUFO HDF5 Flat Field Correction Test Pipeline\n');
fprintf('=============================================\n\n');

% Configuration - Create Config object with proper structure
configData = struct();
configData.output = struct('directory', '/tmp/mufo_test', 'format', 'tiff');
configData.ufo = struct('executable', 'ufo-launch', 'deviceType', 'gpu');
configData.processing = struct('verbose', true, 'precision', 'single');

% File paths - using actual .nxs files
currentDir = pwd;
dataPath = struct(...
    'projections', fullfile(currentDir, 'flyscan_17033-0001.nxs'), ...  % Large projection file (29GB)
    'flats', fullfile(currentDir, 'flyscan_17034-0001.nxs'), ...        % Flat field file (801MB)
    'darks', fullfile(currentDir, 'flyscan_17035-0001.nxs'), ...        % Dark field file (801MB)
    'output', fullfile(configData.output.directory, 'corrected_projections.h5'));

% HDF5 dataset names (actual paths in .nxs files)
h5Dataset = struct(...
    'projections', '/flyscan_17033/scan_data/Image_image', ...
    'flats', '/flyscan_17034/scan_data/Image_image', ...
    'darks', '/flyscan_17035/scan_data/Image_image');

%% Initialize MUFO components
fprintf('Initializing MUFO components...\n');

% Create work directory
workDir = configData.output.directory;
if ~exist(workDir, 'dir')
    mkdir(workDir);
end

% Create Config object and initialize pipeline
config = mufo.Config();
config.set('data.projections.path', dataPath.projections);
config.set('data.flats.path', dataPath.flats);
config.set('data.darks.path', dataPath.darks);
config.set('data.projections.dataset', h5Dataset.projections);
config.set('data.flats.dataset', h5Dataset.flats);
config.set('data.darks.dataset', h5Dataset.darks);
config.set('output.directory', configData.output.directory);

% Initialize pipeline with proper Config object  
pipeline = mufo.Pipeline(config);
loader = mufo.data.DataLoader();

%% Test 1: Load and Inspect Raw Data
fprintf('\n--- Test 1: Loading and Inspecting Raw Data ---\n');

try
    % Load sample projection
    fprintf('Loading sample projection...\n');
    projInfo = h5info(dataPath.projections, h5Dataset.projections);
    projSize = projInfo.Dataspace.Size;
    fprintf('Projection data size: %dx%dx%d\n', projSize);
    
    % Load middle projection for inspection
    midIdx = round(projSize(3)/2);
    sampleProj = h5read(dataPath.projections, h5Dataset.projections, ...
        [1 1 midIdx], [projSize(1) projSize(2) 1]);
    
    % Load flat and dark averages
    fprintf('Loading flat and dark fields...\n');
    flats = h5read(dataPath.flats, h5Dataset.flats);
    darks = h5read(dataPath.darks, h5Dataset.darks);
    
    % Average flats and darks
    flatAvg = mean(flats, 3);
    darkAvg = mean(darks, 3);
    
    % Quick visualization
    figure('Name', 'Raw Data Inspection', 'Position', [100 100 1400 400]);
    
    subplot(1,3,1);
    imagesc(squeeze(sampleProj)); colormap gray; axis image;
    title(sprintf('Sample Projection #%d', midIdx));
    colorbar;
    
    subplot(1,3,2);
    imagesc(flatAvg); colormap gray; axis image;
    title('Flat Field Average');
    colorbar;
    
    subplot(1,3,3);
    imagesc(darkAvg); colormap gray; axis image;
    title('Dark Field Average');
    colorbar;
    
    drawnow;
    
    fprintf('✓ Data loading successful\n');
    
catch ME
    fprintf('✗ Error loading data: %s\n', ME.message);
    rethrow(ME);
end

%% Test 2: UFO Flat Field Correction Pipeline
fprintf('\n--- Test 2: UFO Flat Field Correction ---\n');

try
    % Create temporary files for UFO processing
    tempFiles = struct(...
        'proj', fullfile(workDir, 'temp_proj.tif'), ...
        'flat', fullfile(workDir, 'temp_flat.tif'), ...
        'dark', fullfile(workDir, 'temp_dark.tif'), ...
        'output', fullfile(workDir, 'temp_corrected.tif'));
    
    % Build UFO pipeline
    fprintf('Building UFO pipeline...\n');
    
    % Method 1: Process single slice through UFO
    fprintf('Processing single slice for testing...\n');
    
    % Save test data as TIFF for UFO
    imwrite(uint16(squeeze(sampleProj)), tempFiles.proj);
    imwrite(uint16(flatAvg), tempFiles.flat);
    imwrite(uint16(darkAvg), tempFiles.dark);
    
    % Create UFO command
    chain = mufo.core.UFOChain();
    
    % Read projection
    chain.addTask(mufo.tasks.Read('path', tempFiles.proj));
    
    % Read dark and average (already averaged, but for consistency)
    chain.addTask(mufo.tasks.Read('path', tempFiles.dark));
    
    % Read flat and average
    chain.addTask(mufo.tasks.Read('path', tempFiles.flat));
    
    % Apply flat field correction
    ffcTask = mufo.tasks.FlatFieldCorrect(...
        'absorption_correct', true, ...
        'fix_nan_and_inf', true);
    chain.addTask(ffcTask);
    
    % Write output
    chain.addTask(mufo.tasks.Write('filename', tempFiles.output));
    
    % Execute
    fprintf('Executing UFO flat field correction...\n');
    tic;
    chain.execute();
    elapsed = toc;
    fprintf('✓ UFO processing completed in %.2f seconds\n', elapsed);
    
    % Load corrected slice
    correctedSlice = imread(tempFiles.output);
    
    % Visualize correction result
    figure('Name', 'Flat Field Correction Result', 'Position', [100 100 1200 500]);
    
    subplot(1,3,1);
    imagesc(squeeze(sampleProj)); colormap gray; axis image;
    title('Original Projection');
    caxis([min(sampleProj(:)) max(sampleProj(:))]);
    colorbar;
    
    subplot(1,3,2);
    imagesc(correctedSlice); colormap gray; axis image;
    title('Corrected Projection');
    colorbar;
    
    subplot(1,3,3);
    imagesc(correctedSlice - squeeze(sampleProj)); colormap(gca, 'jet'); axis image;
    title('Difference Map');
    colorbar;
    
    drawnow;
    
catch ME
    fprintf('✗ Error in UFO processing: %s\n', ME.message);
    rethrow(ME);
end

%% Test 3: Full Volume Processing with Progress Monitoring
fprintf('\n--- Test 3: Full Volume Processing ---\n');

try
    % For full volume, we'll process in chunks through UFO
    fprintf('Setting up full volume processing...\n');
    
    % Create output HDF5 file
    if exist(dataPath.output, 'file')
        delete(dataPath.output);
    end
    
    % Process configuration
    chunkSize = 100; % Process 100 projections at a time
    numChunks = ceil(projSize(3) / chunkSize);
    
    % Pre-create HDF5 output file
    h5create(dataPath.output, h5Dataset.projections, projSize, ...
        'Datatype', 'single', ...
        'ChunkSize', [projSize(1) projSize(2) min(chunkSize, projSize(3))], ...
        'Deflate', 4);
    
    % Add metadata
    mufo.data.MetadataManager.writeProcessingInfo(dataPath.output, struct(...
        'process', 'flat_field_correction', ...
        'timestamp', datetime('now'), ...
        'ufo_version', 'latest', ...
        'parameters', ffcTask.parameters));
    
    % Process in chunks
    fprintf('Processing %d chunks of %d projections each...\n', numChunks, chunkSize);
    
    progressBar = waitbar(0, 'Processing projections...', ...
        'Name', 'UFO Flat Field Correction');
    
    totalTime = 0;
    for i = 1:numChunks
        startIdx = (i-1) * chunkSize + 1;
        endIdx = min(i * chunkSize, projSize(3));
        currentChunkSize = endIdx - startIdx + 1;
        
        % Update progress
        waitbar((i-1)/numChunks, progressBar, ...
            sprintf('Processing projections %d-%d of %d', startIdx, endIdx, projSize(3)));
        
        % Load chunk
        projChunk = h5read(dataPath.projections, h5Dataset.projections, ...
            [1 1 startIdx], [projSize(1) projSize(2) currentChunkSize]);
        
        % Process through UFO (simplified for demo - in practice, use UFO's batch processing)
        tic;
        correctedChunk = zeros(size(projChunk), 'single');
        
        % Apply flat field correction formula (simplified version)
        for j = 1:currentChunkSize
            proj = projChunk(:,:,j);
            corrected = (proj - darkAvg) ./ (flatAvg - darkAvg + eps);
            corrected(corrected < 0) = 0;
            corrected(isnan(corrected) | isinf(corrected)) = 0;
            correctedChunk(:,:,j) = corrected;
        end
        
        chunkTime = toc;
        totalTime = totalTime + chunkTime;
        
        % Write chunk to output
        h5write(dataPath.output, h5Dataset.projections, correctedChunk, ...
            [1 1 startIdx], [projSize(1) projSize(2) currentChunkSize]);
        
        % Update estimated time
        avgTimePerChunk = totalTime / i;
        remainingTime = avgTimePerChunk * (numChunks - i);
        waitbar(i/numChunks, progressBar, ...
            sprintf('Processing... ETA: %.1f minutes', remainingTime/60));
    end
    
    close(progressBar);
    fprintf('✓ Full volume processing completed in %.2f minutes\n', totalTime/60);
    
catch ME
    if exist('progressBar', 'var') && isvalid(progressBar)
        close(progressBar);
    end
    fprintf('✗ Error in volume processing: %s\n', ME.message);
    rethrow(ME);
end

%% Test 4: Load and Inspect Corrected Data
fprintf('\n--- Test 4: Loading and Inspecting Corrected Data ---\n');

try
    % Initialize visualization components
    inspector = mufo.vis.Inspector();
    sliceViewer = mufo.vis.SliceViewer();
    
    % Load corrected data info
    correctedInfo = h5info(dataPath.output, h5Dataset.projections);
    fprintf('Corrected data size: %dx%dx%d\n', correctedInfo.Dataspace.Size);
    
    % Load sample slices for inspection
    fprintf('Loading sample slices for inspection...\n');
    
    % Load slices at different positions
    sliceIndices = round(linspace(1, projSize(3), 5));
    sampleSlices = zeros(projSize(1), projSize(2), length(sliceIndices), 'single');
    
    for i = 1:length(sliceIndices)
        sampleSlices(:,:,i) = h5read(dataPath.output, h5Dataset.projections, ...
            [1 1 sliceIndices(i)], [projSize(1) projSize(2) 1]);
    end
    
    % Create comparison figure
    figure('Name', 'Corrected Data Inspection', 'Position', [50 50 1600 900]);
    
    % Show multiple corrected projections
    for i = 1:length(sliceIndices)
        subplot(2, 3, i);
        imagesc(sampleSlices(:,:,i)); colormap gray; axis image;
        title(sprintf('Corrected Projection #%d', sliceIndices(i)));
        colorbar;
    end
    
    % Add histogram comparison
    subplot(2, 3, 6);
    hold on;
    histogram(squeeze(sampleProj), 100, 'Normalization', 'probability', ...
        'DisplayName', 'Original', 'FaceAlpha', 0.5);
    histogram(sampleSlices(:,:,3), 100, 'Normalization', 'probability', ...
        'DisplayName', 'Corrected', 'FaceAlpha', 0.5);
    xlabel('Intensity');
    ylabel('Probability');
    title('Intensity Distribution Comparison');
    legend('Location', 'best');
    grid on;
    
    drawnow;
    
    % Use inspector for detailed view
    fprintf('Opening interactive inspector...\n');
    inspector.inspect(sampleSlices(:,:,3), 'Flat Field Corrected Projection');
    
    fprintf('✓ Data inspection completed\n');
    
catch ME
    fprintf('✗ Error in data inspection: %s\n', ME.message);
    rethrow(ME);
end

%% Test 5: Quality Metrics
fprintf('\n--- Test 5: Computing Quality Metrics ---\n');

try
    % Load original and corrected for comparison
    origProj = h5read(dataPath.projections, h5Dataset.projections, ...
        [1 1 midIdx], [projSize(1) projSize(2) 1]);
    corrProj = h5read(dataPath.output, h5Dataset.projections, ...
        [1 1 midIdx], [projSize(1) projSize(2) 1]);
    
    % Compute metrics
    metrics = struct();
    
    % Signal-to-noise ratio improvement
    noiseRegion = origProj(1:100, 1:100); % Assuming corner is background
    origNoise = std(noiseRegion(:));
    corrNoise = std(corrProj(1:100, 1:100));
    metrics.noiseReduction = (1 - corrNoise/origNoise) * 100;
    
    % Dynamic range
    metrics.origDynamicRange = range(origProj(:));
    metrics.corrDynamicRange = range(corrProj(:));
    
    % Uniformity improvement (coefficient of variation in flat regions)
    flatRegion = corrProj(900:1100, 900:1100); % Center region
    metrics.uniformity = std(flatRegion(:)) / mean(flatRegion(:));
    
    % Display metrics
    fprintf('\nQuality Metrics:\n');
    fprintf('----------------\n');
    fprintf('Noise Reduction: %.1f%%\n', metrics.noiseReduction);
    fprintf('Original Dynamic Range: %.0f\n', metrics.origDynamicRange);
    fprintf('Corrected Dynamic Range: %.0f\n', metrics.corrDynamicRange);
    fprintf('Uniformity (CV): %.4f\n', metrics.uniformity);
    
    % Create metrics visualization
    metricsDisplay = mufo.vis.MetricsDisplay();
    metricsDisplay.showMetrics(metrics, 'Flat Field Correction Quality');
    
    fprintf('\n✓ Quality metrics computed successfully\n');
    
catch ME
    fprintf('✗ Error computing metrics: %s\n', ME.message);
    rethrow(ME);
end

%% Test 6: Memory and Performance Analysis
fprintf('\n--- Test 6: Performance Analysis ---\n');

try
    perfProfiler = mufo.utils.PerformanceProfiler();
    memMonitor = mufo.utils.MemoryMonitor();
    
    % Start monitoring
    memMonitor.start();
    perfProfiler.startTimer('full_pipeline');
    
    % Simulate a complete slice processing
    fprintf('Running performance test on single slice...\n');
    
    perfProfiler.startTimer('data_loading');
    testProj = h5read(dataPath.projections, h5Dataset.projections, ...
        [1 1 1], [projSize(1) projSize(2) 1]);
    perfProfiler.endTimer('data_loading');
    
    perfProfiler.startTimer('preprocessing');
    % Simulate preprocessing
    testCorrected = (testProj - darkAvg) ./ (flatAvg - darkAvg + eps);
    perfProfiler.endTimer('preprocessing');
    
    perfProfiler.startTimer('data_writing');
    tempOutput = fullfile(workDir, 'perf_test.h5');
    h5create(tempOutput, '/data', size(testCorrected), 'Datatype', 'single');
    h5write(tempOutput, '/data', testCorrected);
    perfProfiler.endTimer('data_writing');
    
    perfProfiler.endTimer('full_pipeline');
    
    % Get memory usage
    memStats = memMonitor.stop();
    
    % Display performance results
    fprintf('\nPerformance Summary:\n');
    fprintf('-------------------\n');
    perfProfiler.displayResults();
    
    fprintf('\nMemory Usage:\n');
    fprintf('Peak Memory: %.2f MB\n', memStats.peakMemory / 1e6);
    fprintf('Average Memory: %.2f MB\n', memStats.avgMemory / 1e6);
    
    % Estimate full volume processing time
    sliceTime = perfProfiler.getTime('full_pipeline');
    estimatedTotal = sliceTime * projSize(3);
    fprintf('\nEstimated full volume processing time: %.2f hours\n', estimatedTotal/3600);
    
    fprintf('\n✓ Performance analysis completed\n');
    
catch ME
    fprintf('✗ Error in performance analysis: %s\n', ME.message);
    rethrow(ME);
end

%% Cleanup
fprintf('\n--- Cleanup ---\n');
if exist('tempFiles', 'var')
    fields = fieldnames(tempFiles);
    for i = 1:length(fields)
        if exist(tempFiles.(fields{i}), 'file')
            delete(tempFiles.(fields{i}));
        end
    end
end

fprintf('\n✅ All tests completed successfully!\n');
fprintf('Output saved to: %s\n', dataPath.output);

%% Summary Report
fprintf('\n========== TEST SUMMARY ==========\n');
fprintf('Data Loading: ✓ PASSED\n');
fprintf('UFO Integration: ✓ PASSED\n');
fprintf('HDF5 I/O: ✓ PASSED\n');
fprintf('Visualization: ✓ PASSED\n');
fprintf('Performance: ✓ PASSED\n');
fprintf('==================================\n');