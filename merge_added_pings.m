%merge_added_pings.m
%
%merge missed pings added by fill_gaps.m

clear
close all

import('utilities.readDateTime');

PATH2OUTPUT = 'F:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\GBK_2020_09';
PATH2ADDEDPINGS = 'F:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\GBK_2020_09\ADDED_120kHz_PINGS.mat';
PATH2DETECTIONS = 'F:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\GBK_2020_09\GBK_EK60_DETECTIONS_FILTERED_120kHz_VALIDATED.mat';

load(PATH2ADDEDPINGS);
load(PATH2DETECTIONS);

%load new_pings as table
sz = [length(new_pings) 3];
varTypes = ["string","double","double"];
varNames = ["WavFiles","adjusted_ping_loc","validated"];
new_pings1 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
new_pings1.WavFiles = new_pings(:,1);
new_pings1.adjusted_ping_loc = double(new_pings(:,2));
new_pings1.validated = double(new_pings(:,3));

detections = Filtered_peaks_wav_reviewed(:,{'WavFiles','adjusted_ping_loc','validated'});
unique_wav = unique(Filtered_peaks_wav_reviewed.WavFiles);
merged_detections = [];

for d = 1:length(unique_wav) %start file loop
    file = unique_wav(d);
    dt_start = readDateTime(char((unique_wav(d)))); %start time of file, read in from filename
    new_pings1_sub = new_pings1(new_pings1.WavFiles == file,:);
    detection_subset = detections(detections.WavFiles == file,:);
    merged = [new_pings1_sub;detection_subset];
    merged = sortrows(merged,'adjusted_ping_loc');
    merged_detections = [merged_detections;merged];
end %end file loop

temp = PATH2DETECTIONS(1:end-4);
PATH2SAVE = [temp '_filled.mat'];
save(PATH2SAVE,'merged_detections');