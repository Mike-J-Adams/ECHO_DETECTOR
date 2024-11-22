%combine_validations.m

clear
close all

import('utilities.readDateTime');

%set path to data
PATH2VALIDATIONS = 'F:\BW_ECHO_EXPERIMENT\FCH_2020_09\';

files = '*_VALIDATED.mat';
VALIDATEDList = dir(fullfile(PATH2VALIDATIONS,files));
VALIDATEDFiles = string({VALIDATEDList.name})';
VALIDATED_PINGS = [];

for f = 1:length(VALIDATEDFiles)%start filelist loop
 
    index = VALIDATEDFiles(f);
    temp = split(index,'_');

    PATH2INDEX = char(fullfile(PATH2VALIDATIONS,index));
    disp(PATH2INDEX);    
         %load detections
    load(PATH2INDEX);
   VALIDATED_PINGS = [VALIDATED_PINGS;Filtered_peaks_wav_reviewed];
end

output_name = 'FCH_EK60_VALIDATED_ALL.mat';

save(fullfile(PATH2VALIDATIONS,output_name), "VALIDATED_PINGS")