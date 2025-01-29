%subset_extremes.m

%select and review extreme values

clear
close all

%%%%%%%
PATH2WAV = "D:\BW_ECHO_EXPERIMENT\WAV_FILES\";
WavFileList = dir(fullfile(PATH2WAV,'**\*.wav'));
bandpass_width = 5000;
ping_window = 2048;
%%%%%%%

%load subsetted pings
load("D:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\ALL_SITES\ALL_SITES_PINGS_SPL_filtered.mat")

%only look at PINGS over 5000m
Far_PINGS = ALL_SITES_PINGS_filtered(ALL_SITES_PINGS_filtered.slantRange >= 5000,:);

%%% create empty variable to store bandpass filter object
bandpass_filter = [];
filterFreq = [];

for p = 1:height(Far_PINGS) % loop through pings
    ping_index = Far_PINGS.adjusted_ping_loc(p); % get index of ping peak
    file = Far_PINGS.WavFiles(p); % get wav file
    file_index = WavFileList(find(strcmp(char(WavFileList.name), file )),:);
    filePATH = fullfile(file_index(1).folder,file_index(1).name);
    ainfo = audioinfo(filePATH);
    window_start = ping_index - ping_window/2; % define window to extract ping
    if window_start <= 0
        window_start = 1;
    end
    window_stop = ping_index + ping_window/2;
    if window_stop > ainfo.TotalSamples
        window_stop = ainfo.TotalSamples;
    end
    [x,Fs] = audioread(filePATH,[window_start,window_stop-1]); %read in just ping
    x = x(:,1);
    [M,q] = size(x); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    
    freq = double(Far_PINGS.freq(p)); % get frequency of detected ping
    %define bandpass filter
    LowerStopbandFrequency = freq - bandpass_width-1000;
    LowerPassbandFrequency = freq - bandpass_width;
    UpperPassbandFrequency = freq + bandpass_width;
    UpperStopbandFrequency = freq + bandpass_width+1000;
    
    %bandpass_filter = designfilt(...
    %             'bandpassfir',...
    %             'StopbandFrequency1', LowerStopbandFrequency,...
    %             'PassbandFrequency1', LowerPassbandFrequency,...
    %             'PassbandFrequency2', UpperPassbandFrequency,...
    %             'StopbandFrequency2', UpperStopbandFrequency,...
    %             'StopbandAttenuation1', 60,...
    %             'StopbandAttenuation2', 60,...
    %             'PassbandRipple', 1,...
    %             'DesignMethod', 'kaiserwin',...
    %             'SampleRate', Fs...
    %             ); 
   %Bandpass data
   %x_freq = bandpass.noDelayFilt(bandpass_filter, x1);
   x_freq = bandpass(x,[LowerPassbandFrequency,UpperPassbandFrequency],Fs); %bandpass filter
   
   figure(1)
   subplot(2,1,1)
   plot(x_freq)
       
   figure(1)
   subplot(2,1,2)
   spectrogram(x_freq,64,32,[],Fs,'yaxis')
   colorbar off
   ylim([LowerStopbandFrequency/1000 UpperStopbandFrequency/1000])
      
    
    ui = input('Enter 1 to validated ping, any other key to continue:','s');
    if ui == '1'
        Far_PINGS.validated(p) = 1;
    elseif ui =='0'
        Far_PINGS.validated(p) = 0;
    end  
    pause()
end    

save(fullfile(PATH2WAV,'ALL_SITES_PINGS_SPL_subsetchecked.mat'), "Far_PINGS")





