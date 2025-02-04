%calc_Ping_SPL.m

clear
close all
%%%%%%%
PATH2WAV = "D:\BW_ECHO_EXPERIMENT\WAV_FILES\";
WavFileList = dir(fullfile(PATH2WAV,'**\*.wav'));
bandpass_width = 5000;
ping_window = 2048;
%%%%%%%

%load validated pings
load("D:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\ALL_SITES\ALL_SITES_PINGS.mat")

%load hydrophone values
AMARINPUT = readtable("D:\BW_ECHO_EXPERIMENT\MAPPING_INPUT\AMAR_INPUT.csv");

%HydrophoneSensitivity = -165.42; %C00044
%HydrophoneSensitivity = -165.42; %D000893
%HydrophoneSensitivity = -164.271; %F00099

%%% create empty variable to store bandpass filter object
bandpass_filter = [];
filterFreq = [];

for p = 1:height(ALL_SITES_PINGS) % loop through pings
    ping_index = ALL_SITES_PINGS.adjusted_ping_loc(p); % get index of ping peak
    file = ALL_SITES_PINGS.WavFiles(p); % get wav file
    Site = ALL_SITES_PINGS.SITE(p,:);
    AMAR = AMARINPUT.AMAR(string(AMARINPUT.SITE) == string(Site));
    HSens = AMARINPUT.HydroSens(string(AMARINPUT.SITE) == string(Site));
    Gain = AMARINPUT.Gain(string(AMARINPUT.SITE) == string(Site));
    GainV = AMARINPUT.V(string(AMARINPUT.SITE) == string(Site));
    file_index = WavFileList(find(strcmp(char(WavFileList.name), file )),:);
    filePATH = fullfile(file_index(1).folder,file_index(1).name);
    ainfo = audioinfo(filePATH);
    window_start = ping_index - ping_window; % define window to extract ping
    if window_start <= 0
        window_start = 1;
    end
    window_stop = ping_index + ping_window/2;
    if window_stop > ainfo.TotalSamples
        window_stop = ainfo.TotalSamples;
    end
    [x,Fs] = audioread(filePATH,[window_start,window_stop-1],'native'); %read in just ping in native
    x = double(x(:,1)); %still don't understand this one...
    x1=x/256; %or this one...
    [M,q] = size(x1); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    
    freq = double(ALL_SITES_PINGS.freq(p)); % get frequency of detected ping
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
   x_freq = bandpass(x1,[LowerPassbandFrequency,UpperPassbandFrequency],Fs); %bandpass filter
   trim = 256;
   x_freq_trimmed = x_freq(trim:length(x_freq)-trim); %tring away buffer
   
   figure(1)
   subplot(2,1,1)
   plot(x1)
   subplot(2,1,2)
   plot(x_freq_trimmed)
   
   noise = x_freq_trimmed(1:512); %select out noise
   
   [Start90, Stop90] = utilities.calcEng(x_freq_trimmed, 90); %mark 90% energy 
   hold on
   xline(Start90,'r')
   xline(Stop90,'r')
   hold off
   
    % get min and max magnitudes
    xSub = x_freq_trimmed(Start90:Stop90);
    dur90 = Stop90 - Start90;
    
    avepow = @(x) sum((x.^2),1)./size(x,1); % equation for calculating RMS-based average power for each channel
    pSigWin = avepow(xSub);
    pNoiseWin = avepow(noise);
    ALL_SITES_PINGS.snr_dB(p) = 10*log10(pSigWin./pNoiseWin);
    ALL_SITES_PINGS.snr_adjusted_dB(p) = 10*log10((pSigWin - pNoiseWin)./pNoiseWin);
    
    xMax = max(xSub);
    xMin = min(xSub);
    ppCount = xMax - xMin; %get peak to peak in AMAR counts
    ALL_SITES_PINGS.ppCount(p) = ppCount;
    
   ALL_SITES_PINGS_p = (ppCount*(9/2^24))/(GainV*(10^(HSens/20))); %convert counts to pressure
 
    
    % convert to dB
    ALL_SITES_PINGS.ppSignal(p) = 20*log10(ALL_SITES_PINGS_p); %convert pressure to dB
           
end    

save("D:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\ALL_SITES\ALL_SITES_PINGS_SPL.mat", "ALL_SITES_PINGS")