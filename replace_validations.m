%replace_validations.m

clear
close all


PATH2WAV = "F:\BW_ECHO_EXPERIMENT\";


%load subsetted pings
load("F:\BW_ECHO_EXPERIMENT\ALL_SITES_PINGS_SPL_checked.mat")
load("F:\BW_ECHO_EXPERIMENT\ALL_SITES_PINGS_SPL_subsetchecked.mat")

for i = 1:height(Far_PINGS)
index = find(ALL_SITES_PINGS_filtered_checked.file == Far_PINGS.file(i) & ALL_SITES_PINGS_filtered_checked.adjusted_ping_loc == Far_PINGS.adjusted_ping_loc(i));
ALL_SITES_PINGS_filtered_checked.validated(index) = Far_PINGS.validated(i);
end

ALL_SITES_PINGS_filtered_checked = ALL_SITES_PINGS_filtered_checked(ALL_SITES_PINGS_filtered_checked.validated == 1,:);

%save(fullfile(PATH2WAV,'ALL_SITES_PINGS_SPL_checked.mat'), "ALL_SITES_PINGS_filtered_checked")
%writetable(ALL_SITES_PINGS_filtered_checked,fullfile(PATH2WAV,'ALL_SITES_PINGS_SPL_checked.csv'))