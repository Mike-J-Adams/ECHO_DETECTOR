%review_detections.m
% Review detection and assign 1 or 0 to mask out bad detections

clear
close all

Fs = 256000;
%N=1024;                % sample length
N=Fs*2;                 %N = Fs pull of 2 second of data
N2=N/2;                 % step size = 50% overlap .5 seconds
frequencies =[18000, 38000, 50000, 70000, 120000];
bandpass_width = 2500;  % +- width of bandpass filter

%%% Change to review single frequency, leave empty to review all
freq = 120000;
%%%


PATH2DETECTIONS = 'E:\BW_ECHO_EXPERIMENT\COC_2020_09\OUTPUT';
PATH2DATA = 'E:\BW_ECHO_EXPERIMENT\COC_2020_09\3DaySubset';

files = ['*',int2str(freq),'*.mat'];
DetectionList = dir(fullfile(PATH2DETECTIONS,files));
DetectionFiles = string({DetectionList.name})';

FileList = [];
 
for f = 1:length(DetectionList) %start detection.mat loop
    junk = split(DetectionList(f).name, "_");
    File = junk(5);
    File = strrep(File,'.mat','.wav');
    FileList = [FileList;File];    
end  %end detection.mat loop

FileList = unique(FileList);

for f = 116:length(FileList)%start filelist loop

    wav = FileList(f);
    index = DetectionFiles(f);
    PATH2WAV = char(fullfile(PATH2DATA,wav));
    PATH2INDEX = char(fullfile(PATH2DETECTIONS,index));
    dt_start = readDateTime(PATH2WAV); %start time of file, read in from filename
    [x] = audioread(PATH2WAV); %read in wav file 
    disp(PATH2WAV);
    disp(PATH2INDEX);
    [M,q] = size(x); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    x = detrend(x); %remove mean from audio %%%

    plot_switch1 = 1; %turns test plots on (1) or off (0)
    if plot_switch1 == 1
        figure(1)
        plot(x) %plot entire audio wav with mean removed
        title("Entire audio wav with mean removed")
    end
    
    freq_bins = [freq-bandpass_width freq+bandpass_width]; %create frequency bands using frequencies and width of bandpass
    x_freq = bandpass(x,freq_bins,Fs); %bandpass to isolate each frequency
    Ms = max(abs(x_freq)); %gets maximum from bandpassed audio
    if Ms >= .95
       disp("WARNING: Audio saturated")
    end
    normalx_freq = x_freq/Ms(1); %normalize bandpassed audio %needed?
        
    if plot_switch1 == 1
       figure(2)
       subplot(2,1,1)
       plot(normalx_freq) %plot bandpassed audio wav with mean removed and normalized to max
       title("Bandpassed audio wav")
    end
        
    if plot_switch1 == 1
       figure(2)
       subplot(2,1,2)
       spectrogram(x_freq,16384*2,512,[],Fs,'yaxis')
       colorbar off
       ylim(freq_bins/1000)
    end
    
    %load detections
    load(PATH2INDEX);
        
    freq_detect_SNR.time = freq_detect_SNR.peak_index*dt;
        
    if plot_switch1 == 1
       figure(2)
       subplot(2,1,1)
       hold on
       plot(freq_detect_SNR.peak_index, max(abs(normalx_freq)*1.05), 'r^')
       hold off
    end
    if plot_switch1 == 1
       figure(2)
       subplot(2,1,2)
       hold on
       plot(freq_detect_SNR.time/60,freq/1000, 'r^')
       hold off
    end
        
   start_ping = 0;
   end_ping = 0;
        
   for d = 1:length(freq_detect_SNR.peak_index)
       ping_loc = freq_detect_SNR.peak_index(d);
       ping_window = ping_loc - 5*Fs:ping_loc + 5*Fs-1;
       if ping_window(1) <=0 
          ping_window = 1: 10*Fs;
       end
       if ping_window(end) > length(x_freq)
          ping_window = length(x_freq) - 10*Fs:length(x_freq); 
       end
       start_ping = min(ping_window);
       end_ping = max(ping_window);
       ping = x_freq(ping_window);
       Max_ping = max(abs(ping));
       normal_ping = ping/Max_ping(1);
       [pM,pq] = size(normal_ping); %get size length of audio
       pt = dt*(start_ping:start_ping+pM-1)';
            
       other_pings = freq_detect_SNR.peak_index(freq_detect_SNR.peak_index <= max(ping_window) & freq_detect_SNR.peak_index >= min(ping_window));
            
       figure(4)
       subplot(2,1,1)
       plot(pt,normal_ping)
       hold on
       plot(other_pings*dt,max(normal_ping)*1.06, 'b^')
       plot(ping_loc*dt,max(normal_ping)*1.05, 'r^')
       hold off
       xlim([min(pt) max(pt)])
       subplot(2,1,2)
       spectrogram(normal_ping,16384*2,512,[],Fs,'yaxis')
       colorbar off
       ylim(freq_bins/1000)
           
       pause()
       clf(4)
       clear other_pings ping ping_window
   end %end ping list
end                             %end filelist loop


