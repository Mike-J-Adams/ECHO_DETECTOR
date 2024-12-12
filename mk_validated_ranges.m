%plot_pings.m

clear
close all

import('utilities.readDateTime');

fs = 256000;


%Load validated pings
load("F:\BW_ECHO_EXPERIMENT\COC_2020_09\COC_EK60_VALIDATED_ALL.mat");
validated = VALIDATED_PINGS(VALIDATED_PINGS.validated == 1,:);

%group by minute
validated.adjusted_ping_loc_time = seconds(validated.adjusted_ping_loc/fs);
for i = 1:height(validated)
validated.validation_time(i) = readDateTime(char(validated.WavFiles(i)))+validated.adjusted_ping_loc_time(i);
end
validated.validation_time.Format = 'dd-MMM-uuuu HH:mm:ss.SSSSSS';

validated.validation_time_min = dateshift(validated.validation_time, 'start', 'minute');
figure(1)
subplot(2,1,1)
histogram(validated.validation_time_min,'BinWidth',minutes(1))
xlim([min(validated.validation_time_min),max(validated.validation_time_min)])

Freq_filter = 18000;

validated_freq = validated(str2double(validated.freq) == Freq_filter,:);
subplot(2,1,2)
histogram(validated_freq.validation_time_min,'BinWidth',minutes(1))
xlim([min(validated.validation_time_min),max(validated.validation_time_min)])

% import GPS_min
%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["datetime_1min", "Latitude", "Longitude", "COCranges", "GBKranges", "FCHranges", "COC_SLANT_RANGE", "GBK_SLANT_RANGE", "FCH_SLANT_RANGE"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "datetime_1min", "InputFormat", "yyyy-MM-dd HH:mm:ss");

% Import the data
ALLBIGELOWGPSRANGES = readtable("C:\Users\Adamsmi\Documents\MATLAB\BW_ECHO_EXPERIMENT\ALL_BIGELOW_GPS_RANGES.csv", opts);


% Clear temporary variables
clear opts

%combine pings and ranges
Validated_ranges = join(validated,ALLBIGELOWGPSRANGES,'LeftKeys',18,'RightKeys',1);
save('COC_EK60_VALIDATED_RANGES.mat','Validated_ranges')

Freq_ranges = Validated_ranges(str2double(Validated_ranges.freq) == Freq_filter,:); 
figure(2) 
subplot(3,1,1)
scatter(Freq_ranges.FCHranges, Freq_ranges.FreqRMS90)
title('Outliers included')
xlabel('Horizontal Range (m)')
ylabel('RMS')
%filter out outliers
Freq_ranges_rmoutliers = Freq_ranges(Freq_ranges.FCHranges <=inf,:);
subplot(3,1,2)
scatter(Freq_ranges_rmoutliers.FCHranges, Freq_ranges_rmoutliers.FreqRMS90)
title('Outliers removed')
xlabel('Horizontal Range (m)')
ylabel('RMS')
subplot(3,1,3)
scatter(Freq_ranges_rmoutliers.FCH_SLANT_RANGE, Freq_ranges_rmoutliers.FreqRMS90)
title('Outliers removed')
xlabel('Slant Range (m)')
ylabel('RMS')



