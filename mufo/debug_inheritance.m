% Test inheritance
try
    fprintf('Creating UFOChain directly...\n');
    chain = mufo.core.UFOChain();
    fprintf('UFOChain created successfully\n');
    
    fprintf('Checking if command exists...\n');
    if isprop(chain, 'command')
        fprintf('command property exists\n');
    else
        fprintf('command property does NOT exist\n');
    end
    
    fprintf('Chain properties:\n');
    props = properties(chain);
    for i = 1:length(props)
        fprintf('  %s\n', props{i});
    end
    
catch ME
    fprintf('Error: %s\n', ME.message);
    for i = 1:length(ME.stack)
        fprintf('  %s at line %d\n', ME.stack(i).file, ME.stack(i).line);
    end
end