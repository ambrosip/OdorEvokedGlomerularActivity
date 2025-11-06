function new_file_list = remove_stupid_mac_hidden_files(file_list)

    % wrote by PA with VA and AI help cuz fml

    % Filter out files that start with '._'
    % 1. Get a cell array of just the file names
    file_names = {file_list.name};
    
    % 2. Create a logical index of files that *do not* start with '._'
    % (The tilde '~' negates the result of startsWith)
    files_to_keep_idx = ~startsWith(file_names, '._');
    
    % 3. Use the logical index to keep only the desired entries in the structure array
    new_file_list = file_list(files_to_keep_idx);

end