% Simple debug test
try
    fprintf('Creating Config...\n');
    config = mufo.Config();
    fprintf('Config created successfully\n');
    
    fprintf('Creating Pipeline...\n');
    pipeline = mufo.Pipeline(config);
    fprintf('Pipeline created successfully\n');
    
catch ME
    fprintf('Error: %s\n', ME.message);
    for i = 1:length(ME.stack)
        fprintf('  %s at line %d\n', ME.stack(i).file, ME.stack(i).line);
    end
end