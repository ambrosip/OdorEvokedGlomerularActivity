% plot_z_datasets.m
%
% This script analyzes the uploaded HDF5 file 'analysis_outputs.h5',
% searches for all datasets whose names contain the character 'z',
% reads their data, and plots them in separate figures.

function ai_plot_z_datasets()
    % --- Configuration ---
    filename = '/Volumes/T7 Shield/From Server/20251001/sid237/e1/processed/outputs/batch_analysis_results.hdf5';
    figure_count = 0;
    
    % Check if the file exists before proceeding
    if ~exist(filename, 'file')
        error('HDF5Plotter:FileNotFound', 'Error: The file ''%s'' was not found. Please ensure it is in the current directory.', filename);
    end

    fprintf('Analyzing HDF5 file: %s\n', filename);

    try
        % Get the structure of the HDF5 file
        info = h5info(filename);

        % Start the recursive search from the root of the file
        [figure_count] = traverseAndPlot(info, filename, figure_count);

        if figure_count == 0
            disp('No datasets containing ''z'' were found in the HDF5 file.');
        else
            fprintf('\nFinished plotting. Total of %d dataset(s) containing ''z'' were plotted.\n', figure_count);
        end

    catch ME
        fprintf(2, 'An error occurred during file processing: %s\n', ME.message);
        disp('Check the console for detailed error information.');
    end

end

% --- Helper Function for Recursive Traversal ---
function [figure_count] = traverseAndPlot(info, filename, figure_count)
    
    % 1. Process Datasets in the current group
    for i = 1:length(info.Datasets)
        dataset_name = info.Datasets(i).Name;
        dataset_path = [info.Name, '/', dataset_name];
        
        % Check if the dataset name contains 'z' (case-insensitive)
        if contains(dataset_name, 'z', 'IgnoreCase', true)
            
            fprintf('  -> Found matching dataset: %s\n', dataset_path);
            
            try
                % Read the data
                data = h5read(filename, dataset_path);

                % We primarily want to plot 1D vector data.
                % Check if the data is a 1D vector (i.e., has at most 1 dimension > 1)
                data_size = size(data);
                
                % Determine if it's a plot-able 1D vector (e.g., [N, 1] or [1, N])
                if sum(data_size > 1) <= 1
                    
                    % Flatten the data to ensure correct plotting (e.g., column vector)
                    data = data(:); 
                    
                    % Create a new figure and plot the data
                    figure_count = figure_count + 1;
                    figure(figure_count);
                    plot(data);
                    
                    title(sprintf('Plot of Dataset: %s', dataset_name), 'Interpreter', 'none');
                    xlabel('Index');
                    ylabel('Value');
                    grid on;
                    
                    fprintf('     (Plotted successfully in Figure %d)\n', figure_count);

                elseif sum(data_size > 1) == 2
                    % If it is 2D, attempt a surface plot or image/heatmap
                    fprintf('     (2D data detected. Attempting to plot as a heatmap.)\n');
                    figure_count = figure_count + 1;
                    figure(figure_count);
                    imagesc(data);
                    colorbar;
                    axis tight;
                    title(sprintf('Heatmap of Dataset: %s', dataset_name), 'Interpreter', 'none');
                    xlabel('Column Index');
                    ylabel('Row Index');

                else
                    % Data is not a simple 1D or 2D array, skip plotting
                    fprintf('     (Skipped plotting: Data size is %s, not a simple 1D or 2D array.)\n', mat2str(data_size));
                end

            catch ME_READ
                fprintf(2, '     (Error reading dataset %s: %s)\n', dataset_path, ME_READ.message);
            end
        end
    end
    
    % 2. Process Groups (Recurse into subgroups)
    for j = 1:length(info.Groups)
        fprintf('Entering Group: %s\n', info.Groups(j).Name);
        % Recursive call to continue traversal
        [figure_count] = traverseAndPlot(info.Groups(j), filename, figure_count);
    end
end

% % Execute the main function when the script runs
% ai_plot_z_datasets();