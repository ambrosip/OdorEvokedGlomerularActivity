%{
DOCUMENTATION
written by VA on Apr/2026

GOAL:
    Calculate and save movies of dF/F, z-scores of F (over time),
    and z-score of dF/F (over time and over time and space).

ALERT: 
    1) Odor presentation window is marked with a red circle in the video.
    2) We average across acquisitions before computing dF/F.
    3) Don't save the workspace (it's too big).
    4) Saves to processed > matlab > movies

DEPENDS on:
    mat file from preProcessing
%}

%% USER INPUT

% Experiment folder
expFolder = expDir;

% Do you want to save movies as tiff too?
saveAsTiff = 0;

% Define the colors
max_df_color = [103 0 31] / 255;
min_df_color = [5 48 97] / 255;

% Range of z-score to include in movies
plotRange = [-5, 5];

% How many seconds in the rolling window
secondsRollingWindow = .5;

% Length in seconds for the baseline
% secondsBaselineDuration = 1.0;


%% Creates Diverging Colormap

% Creates 512 color steps between the two colors (512 is arbitrary)
% The number just need to be high enough for the gradient to be smooth
colorpoints = linspace(0.0, 1.0, 512);

% This function is from the following website:
% https://www.kennethmoreland.com/color-maps/
divergingGradient = divergingMap(colorpoints, min_df_color, max_df_color);


%% Compute Timing Info

frameOdorOnset = round(odor_onset_s * frame_rate_hz);
frameOdorOffset = round(odor_offset_s * frame_rate_hz);
frameOdorDuration = frameOdorOffset - frameOdorOnset + 1;

framesRollingWindow = round(secondsRollingWindow * frame_rate_hz);

% ASSUMPTION: 1) Baseline has the same size as rolling window
%             2) Baseline ends just before odor presentation

frameBaselineStart = frameOdorOnset - framesRollingWindow;
frameOdorOnsetFromBaseline = framesRollingWindow + 1;
frameOdorOffsetFromBaseline = framesRollingWindow + frameOdorDuration;


%% Get mcor File Paths and saveFolder

expDir = expFolder;

if exist('patchwarp','var')
    if patchwarp == 1
        getFileDirs_patchwarp
    else
        getFileDirs
    end
else
    getFileDirs
end

getImgDirs

todayStr = string(datetime('Today'), 'yyyy-MM-dd');
analysisDate = todayStr;
saveFolder = fullfile(expFolder, ...
    'processed', 'matlab', 'movies');

if ~isfolder(saveFolder)
    mkdir(saveFolder);
end

% Assumes the default folder structure for a experiment
% mcorFolder = fullfile(expFolder, 'processed', 'mcor');
mcorFolder = mcorImgDir;

% Store file paths
% mcorFilePaths = dir(fullfile(mcorFolder, "*_mcor.tif"));
mcorFilePaths = dir(fullfile(mcorFolder, "*.tif"));
if isempty(mcorFilePaths), error("mcor folder is empty."); end

% Make sure they are in the right order
[~, ind] = sort({mcorFilePaths.name});
mcorFilePaths = mcorFilePaths(ind);


%% Compute Average Signals

fprintf("[INFO] Loading files and computing averages...\n")

% Find the sizes of the movies
% We will not load anything before the baseline
filename = mcorFilePaths(1).name;
filepath = fullfile(mcorFolder, filename);
infoFile = imfinfo(filepath);

framesAfterBaselineStart = length(infoFile) - frameBaselineStart + 1;
movieSize = [infoFile(1).Height ...
             infoFile(1).Width ...
             framesAfterBaselineStart];

% There is one movie per program and per odor
% The index below tracks how many movies we started so far
iMovie = 0;

% Iterate through all the possibilities
for program_name = string(fieldnames(s_olfactometer))'
    program = s_olfactometer.(program_name);

    % Get program number (number after 'program_')
    programSplit = split(program_name, '_');
    programNumber = str2double(programSplit(2));

    % Skip programs that don't have summaries
    if ~isfield(program, 'summary_by_trial')
        message = ['[WARNING] Skipping program %d because ' ...
            'it was marked as ''ignore''.\n'];
        fprintf(message, programNumber)
        continue
    end

    % Get program type
    programType = s_olfactometer.(program_name).type;

    % Get odor and outcome info
    summaryByTrial = program.summary_by_trial;
    programOdors = unique(summaryByTrial.odor);

    for odor = programOdors'
        fprintf("[INFO]   Processing program %s odor %d\n", ...
            programType, odor);

        % One more movie is being added
        iMovie = iMovie + 1;

        % Add movie to list of movies
        movies(iMovie).program = program_name;
        movies(iMovie).type = programType;
        movies(iMovie).odor = odor;
        movies(iMovie).acquisitions = 0;

        % Get table with data only for the current odor
        odorRows = summaryByTrial.odor == odor;
        odorTable = summaryByTrial(odorRows, :);

        % Initialize counts of how many images where used in the average
        movies(iMovie).hit.total = 0;
        movies(iMovie).miss.total = 0;
        movies(iMovie).false.total = 0;
        movies(iMovie).na.total = 0;

        % reverseStr = '';
        % msg = sprintf('[INFO] Acquisitions processed: %3.1f%%\n', 0);
        % fprintf([reverseStr, msg]);
        
        for row = 1:height(odorTable)
            % Process acquisition
            acqIdx = odorTable{row, 'acqIdx'};

            % Skip trials without acquisition
            if isnan(acqIdx) || acqIdx == 0
                continue
            end

            filename = mcorFilePaths(acqIdx).name;
            fileDir = mcorFilePaths(acqIdx).folder;
            filepath = fullfile(fileDir, filename);

            % Skips computation and prints warning if file not found
            if ~isfile(filepath)
                fprintf('WARNING: File %s not found.\n', filename)
                continue
            end

            % Load relevant part of the mcor movie
            movie = single(read_file(filepath, frameBaselineStart)) + 32768;            

            % outcome is the first word of the outcome name
            % Possibilities are "hit", "miss", "false", and "na".
            outcome = split(odorTable{row, 'outcome'});
            outcome = outcome(1);

            % One more acquisition for this outcome
            movies(iMovie).acquisitions = ...
                movies(iMovie).acquisitions + 1;
            movies(iMovie).(outcome).total = ...
                movies(iMovie).(outcome).total + 1;

            n = movies(iMovie).(outcome).total;

            % Only allocate movie if there is an acquisition
            if n == 1
                movies(iMovie).(outcome).avgMovie = zeros(movieSize);
            end

            % Averaging across acquisitions
            % We use a simple trick to compute the average on the fly:
            % https://stackoverflow.com/a/23493727

            movies(iMovie).(outcome).avgMovie = ...
                movies(iMovie).(outcome).avgMovie * (n - 1) / n ...
                + movie / n;

            % Calculate and display progress
            % percentDone = 100 * row / height(odorTable);
            % 
            % msg = sprintf( ...
            %     '[INFO] Acquisitions processed: %3.1f%%\n', percentDone);
            % fprintf([reverseStr, msg]);
            % 
            % reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
        
    end
end

%% Compute Z-Scores for the dF/F of the Averages

fprintf("[INFO] Computing z-score of dF/F...\n")

nMovies = length(movies);

for iMovie = 1:nMovies
    for outcome = ["hit" "miss" "false" "na"]
        % Movie for that outcome exists if there is at least one
        % acquisition with that outcome.
        if movies(iMovie).(outcome).total > 0
            [ movies(iMovie).(outcome).FZScore ...
            , movies(iMovie).(outcome).dFF ...
            , movies(iMovie).(outcome).dFFZScoreTime ...
            , movies(iMovie).(outcome).dFFZScoreAll ...
            ] = utils.movingZScoreAll( ...
                movies(iMovie).(outcome).avgMovie, framesRollingWindow);
        end
    end
end

%% Save movies

fprintf("[INFO] Saving %d movies:\n", nMovies);

% Make all imshow figures borderless
iptsetpref('ImshowBorder', 'tight');

for iMovie = 1:nMovies
    for outcome = ["hit" "miss" "false" "na"]
        for movieType = ["FZScore", "dFF", "dFFZScoreTime", "dFFZScoreAll"]
        
            % Don't create a movie if there is no data
            if movies(iMovie).(outcome).total == 0
                continue
            end
    
            fmt = "%s_to_%s_zScoreMovie_%s_%s_odor_%d_outcome_%s_%s.avi";
            videoName = sprintf( fmt ...
                               , firstAcqName(1:22) ...
                               , lastAcqName(1:22) ...
                               , movies(iMovie).program ...
                               , movies(iMovie).type ...
                               , movies(iMovie).odor ...
                               , outcome ...
                               , movieType);

            videoPath = fullfile(saveFolder, videoName);
            fprintf("[INFO]   Movie %d -> %s\n", iMovie, videoName);
    
            video = VideoWriter(videoPath, "MPEG-4");
            video.FrameRate = frame_rate_hz;
            open(video);
    
            movie = movies(iMovie).(outcome).(movieType);
    
            % Encoding works better if width and height are multiples of 16
            w = 16 * ceil(size(movie, 1) / 16);
            h = 16 * ceil(size(movie, 2) / 16);
    
            fig = figure('Position', [50 50 w h], 'Visible', 'off');
    
            for iFrame = 1:size(movies(iMovie).(outcome).avgMovie, 3)
                imshow(movie(:, :, iFrame), plotRange, ...
                    'InitialMagnification', 'fit');
                colormap(divergingGradient);
    
                if iFrame > frameOdorOnsetFromBaseline && ...
                        iFrame < frameOdorOffsetFromBaseline
    
                    % Draw a circle if inside odor presentation window
                    rectangle('Position', [20 20 10 10], ...
                              'Curvature', [1, 1], ...
                              'FaceColor', 'r');
                end
    
                frame = getframe(fig);
                writeVideo(video, frame);
            end
    
            close(fig);
            close(video);
        end
    end
end

fprintf("[INFO] Finished saving!\n")


%% Save movies as TIFF

if saveAsTiff == 1

    fprintf("[INFO] Saving %d movies:\n", nMovies);
        
    for iMovie = 1:nMovies
        for outcome = ["hit" "miss" "false" "na"]
            % Don't create a movie if there is no data
            if movies(iMovie).(outcome).total == 0
                continue
            end
    
            tiffName = strcat(firstAcqName(1:22), "_to_", ...
                lastAcqName(1:22), "_", ...
                sprintf(...              
                "zScoreMovie_%s_%s_odor_%d_outcome_%s.tif", ...
                movies(iMovie).program, ...
                movies(iMovie).type, ...
                movies(iMovie).odor, ...
                outcome));
            tiffPath = fullfile(saveFolder, tiffName);
            fprintf("[INFO]   Movie %d -> %s\n", iMovie, tiffName);
    
            movie = movies(iMovie).(outcome).avgMovie;
            t = Tiff(tiffPath, "w");
    
            for iFrame = 1:size(movies(iMovie).(outcome).avgMovie, 3)
                tagstruct.ImageLength = size(movie, 1);
                tagstruct.ImageWidth = size(movie, 2);
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
                tagstruct.BitsPerSample = 32;
    
                t.setTag(tagstruct);
                t.write(single(movie(:, :, iFrame)));
    
                if iFrame < size(movie, 3)
                    t.writeDirectory();
                end
            end
    
            t.close();
        end
    end
    
    fprintf("[INFO] Finished saving!\n")
end