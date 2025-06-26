% Explore scan_data group
try
    fprintf('Exploring scan_data group in flyscan_17033-0001.nxs...\n');
    
    % Check scan_data group
    scan_info = h5info('flyscan_17033-0001.nxs', '/flyscan_17033/scan_data');
    
    fprintf('Datasets in scan_data:\n');
    for i = 1:length(scan_info.Datasets)
        dataset = scan_info.Datasets(i);
        fprintf('  %s: size %s, datatype %s\n', dataset.Name, ...
                mat2str(dataset.Dataspace.Size), dataset.Datatype.Class);
    end
    
    fprintf('\nGroups in scan_data:\n');
    for i = 1:length(scan_info.Groups)
        group = scan_info.Groups(i);
        fprintf('  %s/\n', group.Name);
    end
    
    % Look for detector data
    if ~isempty(scan_info.Groups)
        for i = 1:length(scan_info.Groups)
            group = scan_info.Groups(i);
            fprintf('\nExploring %s:\n', group.Name);
            
            if ~isempty(group.Datasets)
                for j = 1:min(10, length(group.Datasets))
                    dataset = group.Datasets(j);
                    fprintf('  %s: size %s\n', dataset.Name, mat2str(dataset.Dataspace.Size));
                end
            end
        end
    end
    
catch ME
    fprintf('Error exploring scan_data: %s\n', ME.message);
end