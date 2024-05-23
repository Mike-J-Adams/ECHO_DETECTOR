%combine_detections.m
%combine detections and filter by SNR and duration

clear
close all

import('utilities.readDateTime');

%set path to data
PATH2DETECTIONS = 'D:\BW_ECHO_EXPERIMENT\COC_2020_09\ECHO_DETECTOR_OUTPUT\OUTPUT';

%set thresholds
SNR_THRESHOLD = 5;
Ping_Duration = [0.001 0.05];

files = '*.mat';
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
    peaks.freq = repmat(temp(4),height(peaks),1);
    peaks.file = repmat(temp(5),height(peaks),1);
    PEAKS=[PEAKS;peaks];
end

Filtered_PEAKS = PEAKS(PEAKS.FreqSNR2>=SNR_THRESHOLD & PEAKS.FreqDUR90>=Ping_Duration(1) & PEAKS.FreqDUR90<=Ping_Duration(2) ,:);
test = Filtered_PEAKS(Filtered_PEAKS.freq == '18000',:);

output_name = 'COC_EK60_DETECTIONS.mat';

save(output_name, "Filtered_PEAKS")