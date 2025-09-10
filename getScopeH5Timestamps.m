%% Get trial start and odor presentation timestamps from h5 file

% get relevant data
total_data_points = h5info(h5_file_dir, trial_start_label_h5).Dataspace.Size;
samplerate = h5info(h5_file_dir).Attributes.Value;
trial_start_TTL = h5read(h5_file_dir, trial_start_label_h5);
odor_TTL = h5read(h5_file_dir, odor_dur_h5);
[~,fileName_h5,~] = fileparts(h5_file_dir);

% x axis in minutes
x_data_points_h5 = 1:total_data_points;
x_minutes_h5 = x_data_points_h5/samplerate/60;

% find onset of TTL pulses
[trial_pks,trial_locs]=findpeaks(diff(trial_start_TTL),'MinPeakHeight',2,'MinPeakDistance',samplerate/2);
    % added MinPeakDistance of 500 ms to deal with problematic file where
    % code found 2 trial_locs right next to each other
[odor_pks,odor_locs]=findpeaks(diff(odor_TTL),'MinPeakHeight',2,'MinPeakDistance',samplerate/10);
    % added MinPeakDistance of 100 ms to deal with problematic file where
    % code found 2 odor_locs right next to each other

% find offset of odor TTL pulse
[odor_end_pks,odor_end_locs]=findpeaks(-diff(odor_TTL),'MinPeakHeight',2);

% adjust timing of locs (to account for diff function used to find peaks)
trial_locs = trial_locs + 1;
odor_locs = odor_locs + 1;
odor_end_locs = odor_end_locs + 1;

% convert locs from data points to minutes
trial_locs = trial_locs/samplerate/60;
odor_locs = odor_locs/samplerate/60;
odor_end_locs = odor_end_locs/samplerate/60;

disp('got h5 data')