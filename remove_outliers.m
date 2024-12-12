%remove_outliers.m

%review and remove outlier validations

clear
close all

%load pings
load("F:\BW_ECHO_EXPERIMENT\ALL_SITES_PINGS_SPL.mat")

%filter validations
% range
ALL_SITES_PINGS_filtered = ALL_SITES_PINGS(ALL_SITES_PINGS.slantRange <= 8000,:); 
ALL_SITES_PINGS_filtered = ALL_SITES_PINGS_filtered(ALL_SITES_PINGS_filtered.snr_dB > 6,:);
ALL_SITES_PINGS_filtered = ALL_SITES_PINGS_filtered(string(ALL_SITES_PINGS_filtered.freq) ~= "50000",:);

%save('ALL_SITES_PINGS_SPL_filtered.mat','ALL_SITES_PINGS_filtered')

figure(1)
scatter(ALL_SITES_PINGS_filtered.snr_dB,ALL_SITES_PINGS_filtered.ppSignal)
title("SNR (dB) vs P-P Signal (dB)")


%break down by site
SITES = unique(string(ALL_SITES_PINGS_filtered.SITE));

figure(2)
for i = 1:length(SITES)
temp = ALL_SITES_PINGS_filtered(string(ALL_SITES_PINGS_filtered.SITE) == SITES(i),:); 
scatter(temp.slantRange,temp.ppSignal)
hold on
end
title('Slant Range')
xlabel('Slant Range (m)')
ylabel('SPL (dB)')
legend(SITES)
hold off


%breakdown by frequency
freq = unique(string(ALL_SITES_PINGS_filtered.freq));

figure(3)
for f = 1:length(freq)
temp = ALL_SITES_PINGS_filtered(string(ALL_SITES_PINGS_filtered.freq) == freq(f),:); 
scatter(temp.slantRange,temp.ppSignal)
hold on
end
title('Frequency')
xlabel('Horizontal Range (m)')
ylabel('SPL (dB)')
legend(freq)
hold off

%Breakdown range v time
figure(4)
for f = 1:length(freq)
temp = ALL_SITES_PINGS_filtered(string(ALL_SITES_PINGS_filtered.freq) == freq(f),:); 
scatter(temp.validation_time_min,temp.slantRange)
hold on
end
title('Frequency')
xlabel('Datetime')
ylabel('Slant Range (m)')
legend(freq)
hold off

%plot cloverleafs

for i = 1:length(SITES)
    temp = ALL_SITES_PINGS_filtered(string(ALL_SITES_PINGS_filtered.SITE) == SITES(i),:); 
    figure(4+i)
    for f = 1:length(freq)
        temp2 = temp(string(temp.freq) == freq(f),:);
        scatter(temp2.Longitude,temp2.Latitude)
        hold on
    end
end