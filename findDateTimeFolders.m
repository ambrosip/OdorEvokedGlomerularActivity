% written by AI with PA's prompt: write matlab code that finds all
% downstream folders from a given directory that have the following
% structure: yyyy_mm_dd-hh_mm_ss 

function matchingFolders = findDateTimeFolders(startDir)
    % findDateTimeFolders finds all subfolders with name format yyyy_mm_dd-hh_mm_ss
    % starting from a given directory.
    %
    % Input:
    %   startDir - The root directory path (e.g., '/home/user/data' or 'C:\Data')
    %
    % Output:
    %   matchingFolders - A cell array of full paths to the matching folders

    if nargin < 1
        startDir = pwd; % Use current directory if none specified
    end

    % Use the '**' wildcard with dir to list all files and folders recursively
    % The trailing '*' is important for dir to list all items within the
    % subdirectories matched by '**'
    all_items = dir(fullfile(startDir, '**', '*'));

    % Filter for directories only and exclude '.' and '..'
    is_dir = [all_items.isdir];
    all_dirs = all_items(is_dir);
    
    % Filter out the '.' and '..' entries which are always present
    valid_dirs_idx = ~strcmp({all_dirs.name}, '.') & ~strcmp({all_dirs.name}, '..');
    valid_dirs = all_dirs(valid_dirs_idx);

    % Define the regular expression for the folder name format: yyyy_mm_dd-hh_mm_ss
    % pattern components:
    % \d{4} - four digits (yyyy)
    % _     - underscore literal
    % \d{2} - two digits (mm, dd, hh, mm, ss)
    % -     - hyphen literal
    regex_pattern = '^\d{4}_\d{2}_\d{2}-\d{2}_\d{2}_\d{2}$';

    % Use regexp to find names that match the pattern
    % 'match' keyword returns the matching substrings in a cell array
    dir_names = {valid_dirs.name};
    matches = regexp(dir_names, regex_pattern, 'match', 'once');

    % Filter out empty cells (non-matches)
    is_match = ~cellfun('isempty', matches);
    matched_dirs = valid_dirs(is_match);

    % Get the full paths of the matching directories
    matchingFolders = cell(1, numel(matched_dirs));
    for i = 1:numel(matched_dirs)
        matchingFolders{i} = fullfile(matched_dirs(i).folder, matched_dirs(i).name);
    end
end

% Example Usage:
% Assuming a folder structure exists within the current working directory.
% start_directory = pwd; 
% matched_list = findDateTimeFolders(start_directory);
% 
% disp('Found matching folders:');
% disp(matched_list');