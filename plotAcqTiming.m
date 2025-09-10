%% FIG 2 - SCOPE H5 Time-series and peaks

% the timing of the 1st acq is 0, so to align acqs to trials, I need to
% discount the time interval between the start of the loop and the start of
% the first trial
adjusted_rawImgStarts_min = rawImgStarts_min + trial_locs(1);

fig2 = figure('name', strcat(fileName_h5, '_', analysisDate, ' - scope events'));
plot(x_minutes_h5,trial_start_TTL, 'Color','k')
hold on;
plot(trial_locs,trial_pks,'o','Color','k')
plot(x_minutes_h5,odor_TTL, 'Color','m')
% plot(odor_locs,trial_pks,'o','Color','m')
plot(odor_locs,odor_pks,'o','Color','m')
plot(odor_end_locs,odor_end_pks,'o','Color','b')
% show acq start times
xline(adjusted_rawImgStarts_min)
xlabel('Time (min)')
ylabel('TTLs (trial start, odor) and Events (acqs)','Interpreter','none')
hold off
disp('plot fig2 complete')
