%Plot_all_pings.m

clear
close all

%load pings
load("F:\BW_ECHO_EXPERIMENT\ALL_SITES_PINGS_SPL_checked.mat");

pings = ALL_SITES_PINGS_filtered_checked;

figure(1)
scatter(pings.snr_dB,pings.ppSignal)
title("SNR (dB) vs P-P Signal (dB)")


%break down by site
SITES = unique(string(pings.SITE));

figure(2)
for i = 1:length(SITES)
temp = pings(string(pings.SITE) == SITES(i),:); 
scatter(temp.slantRange,temp.ppSignal)
hold on
end
title('Slant Range')
xlabel('Slant Range (m)')
ylabel('SPL (dB)')
legend(SITES)
hold off



%breakdown by frequency
freq = unique(string(pings.freq));



figure(3)
for f = 1:length(freq)
temp = pings(string(pings.freq) == freq(f),:); 
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
temp = pings(string(pings.freq) == freq(f),:); 
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
    temp = pings(string(pings.SITE) == SITES(i),:); 
    figure(4+i)
    for f = 1:length(freq)
        temp2 = temp(string(temp.freq) == freq(f),:);
        scatter(temp2.Longitude,temp2.Latitude)
        hold on
    end
end


