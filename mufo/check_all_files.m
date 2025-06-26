% Check all three .nxs files
files = {'flyscan_17033-0001.nxs', 'flyscan_17034-0001.nxs', 'flyscan_17035-0001.nxs'};

for f = 1:length(files)
    fprintf('\n=== %s ===\n', files{f});
    try
        % Get top level structure
        info = h5info(files{f});
        root_group = info.Groups(1).Name(2:end); % Remove leading /
        
        % Check scan_data
        scan_path = ['/' root_group '/scan_data'];
        scan_info = h5info(files{f}, scan_path);
        
        fprintf('Datasets in scan_data:\n');
        for i = 1:length(scan_info.Datasets)
            dataset = scan_info.Datasets(i);
            fprintf('  %s: size %s\n', dataset.Name, mat2str(dataset.Dataspace.Size));
        end
        
        % Check for Image_image dataset specifically
        image_path = [scan_path '/Image_image'];
        try
            image_info = h5info(files{f}, image_path);
            fprintf('Image data path: %s (size: %s)\n', image_path, mat2str(image_info.Dataspace.Size));
        catch
            fprintf('No Image_image dataset found\n');
        end
        
    catch ME
        fprintf('Error reading %s: %s\n', files{f}, ME.message);
    end
end