%combine_filled_detections.m
%combine filled detections and filter by SNR and duration

clear
close all

import('utilities.readDateTime');

%set path to data
PATH2DETECTIONS = 'F:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\GBK_2020_09';

files = '*filled.mat';
DetectionList = dir(fullfile(PATH2DETECTIONS,files));
DetectionFiles = string({DetectionList.name})';
PEAKS = [];

for f = 1:length(DetectionFiles)%start filelist loop
 
    index = DetectionFiles(f);
    temp = split(index,'_');

    PATH2INDEX = char(fullfile(PATH2DETECTIONS,index));
    disp(PATH2INDEX);    
         %load detections
    load(PATH2INDEX);
    PEAKS=[PEAKS;merged_detections];
end

%Filtered_PEAKS = PEAKS(PEAKS.FreqSNR2>=SNR_THRESHOLD & PEAKS.FreqDUR90>=Ping_Duration(1) & PEAKS.FreqDUR90<=Ping_Duration(2) ,:);
%test = Filtered_PEAKS(Filtered_PEAKS.freq == '18000',:);
%test2 = PEAKS(PEAKS.freq == '18000',:);
output_name = 'GBK_EK60_DETECTIONS_ALL.mat';

save(output_name, "PEAKS")