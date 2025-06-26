% Explore .nxs file structure
try
    fprintf('Exploring flyscan_17033-0001.nxs structure...\n');
    info = h5info('flyscan_17033-0001.nxs');
    
    fprintf('Root level groups:\n');
    for i = 1:length(info.Groups)
        fprintf('  /%s\n', info.Groups(i).Name);
    end
    
    % Look for common dataset paths
    common_paths = {'/entry/data/data', '/entry/instrument/detector/data', ...
                   '/entry/instrument/detector/image_key', '/entry/sample/projections'};
    
    fprintf('\nChecking common dataset paths:\n');
    for i = 1:length(common_paths)
        try
            dataset_info = h5info('flyscan_17033-0001.nxs', common_paths{i});
            fprintf('  %s: EXISTS - Size: %s\n', common_paths{i}, mat2str(dataset_info.Dataspace.Size));
        catch
            fprintf('  %s: NOT FOUND\n', common_paths{i});
        end
    end
    
    % Recursively explore the first group
    if ~isempty(info.Groups)
        first_group = info.Groups(1);
        fprintf('\nExploring %s:\n', first_group.Name);
        if ~isempty(first_group.Groups)
            for i = 1:min(5, length(first_group.Groups))
                fprintf('  %s/\n', first_group.Groups(i).Name);
            end
        end
        if ~isempty(first_group.Datasets)
            for i = 1:min(5, length(first_group.Datasets))
                fprintf('  %s (size: %s)\n', first_group.Datasets(i).Name, mat2str(first_group.Datasets(i).Dataspace.Size));
            end
        end
    end
    
catch ME
    fprintf('Error exploring file: %s\n', ME.message);
end