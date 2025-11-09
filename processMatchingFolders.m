% written by AI with PA's prompt: now write another function that will
% iterate through all directories found in matchingFolders and apply the
% function "analyzeBehavior" 

function processMatchingFolders(matchingFolders,saveDir)
    % processMatchingFolders iterates through a cell array of directory paths
    % and applies the 'analyzeBehavior' function to each path.

    if isempty(matchingFolders)
        disp('No matching folders found to process.');
        return;
    end

    disp(['Starting analysis for ' num2str(numel(matchingFolders)) ' folders...']);

    % Iterate through each folder path provided in the input cell array
    for i = 1:numel(matchingFolders)
        current_folder = matchingFolders{i};
        
        fprintf('--- Processing folder %d/%d: %s ---\n', i, numel(matchingFolders), current_folder);
        
        % Ensure the target function 'analyzeBehavior' exists
        if exist('analyzeBehavior', 'file')
            try
                % Call the external function with the current folder path
                analyzeBehavior(current_folder,saveDir);
                disp('Analysis complete for this folder.');
            catch ME
                % Catch potential errors during the analysis of a specific folder
                warning('Error processing folder: %s\nError message: %s', current_folder, ME.message);
            end
        else
            error('The function "analyzeBehavior" does not exist or is not in the MATLAB path.');
        end
    end

    disp('All specified folders have been processed.');
end
