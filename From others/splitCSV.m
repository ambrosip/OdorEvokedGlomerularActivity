% Goal: split csv file into two halves
% Why? split arduino treadmill data between e1 and e2
% written with Gemini 

% Define the name of your original CSV file
filename = wheelDataDir; 

% Read the entire CSV file into a table
T = readtable(filename);

% Determine the number of rows in the table
numRows = size(T, 1);

% Calculate the splitting index (using floor for an integer split point)
splitIndex = floor(numRows / 2);

% Split the table into two halves
table1 = T(1:splitIndex, :);
table2 = T(splitIndex+1:end, :);

% Define names for the new CSV files
filename1 = 'first_half.csv';
filename2 = 'second_half.csv';

% Write the two halves to separate CSV files
fulldirectory1 = fullfile(saveDir, strcat(wheelDataName, '_', filename1));
fulldirectory2 = fullfile(saveDir, strcat(wheelDataName, '_', filename2));
writetable(table1, fulldirectory1, 'WriteMode', 'overwrite');
writetable(table2, fulldirectory2, 'WriteMode', 'overwrite');

disp(['Successfully split ' filename ' into ' fulldirectory1 ' and ' fulldirectory2 '.']);