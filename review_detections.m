%review_detections.m
% Review detection and assign 1 or 0 to mask out bad detections

clear
close all

import('utilities.readDateTime');


plot_switch1 = 1; %turns test plots on (1) or off (0)
Fs = 256000;
%N=1024;                % sample length
N=Fs*2;                 %N = Fs pull of 2 second of data
N2=N/2;                 % step size = 50% overlap .5 seconds
frequencies =[18000, 38000, 50000, 70000, 120000];
bandpass_width = 5000;  % +- width of bandpass filter

SNR_THRESHOLD = 6;
Ping_Duration = [0.001 0.0025];

%%% Change to review single frequency, leave empty to review all
freq = 18000;
%%%

PATH2OUTPUT = 'D:\BW_ECHO_EXPERIMENT\COC_2020_09';
output_name = 'COC_EK60_DETECTIONS_FILTERED_VALIDATED.mat';
PATH2DETECTIONS = 'D:\BW_ECHO_EXPERIMENT\COC_2020_09\COC_EK60_DETECTIONS_FILTERED_VALIDATED.mat';
PATH2DATA = 'D:\BW_ECHO_EXPERIMENT\COC_2020_09\3DaySubset';
load(PATH2DETECTIONS);

if exist('Filtered_PEAKS')
    PEAKS = Filtered_PEAKS;
end

%%% create empty variable to store bandpass filter object
bandpass_filter = [];
filterFreq = [];

Filtered_PEAKS = PEAKS(PEAKS.FreqSNR2>=SNR_THRESHOLD & PEAKS.FreqDUR90>=Ping_Duration(1) & PEAKS.FreqDUR90<=Ping_Duration(2) & PEAKS.freq == num2str(freq) ,:);
Filtered_PEAKS.WavFiles = strrep(Filtered_PEAKS.file,'.mat','.wav');
unique_wav = unique(Filtered_PEAKS.WavFiles);
if ~ismember('validated', Filtered_PEAKS.Properties.VariableNames)
    Filtered_PEAKS.validated = zeros(size(Filtered_PEAKS,1),1);
end
if ~ismember('reviewed', Filtered_PEAKS.Properties.VariableNames)
    Filtered_PEAKS.reviewed = zeros(size(Filtered_PEAKS,1),1);
else
    Filtered_PEAKS = Filtered_PEAKS(Filtered_PEAKS.reviewed == 0,:);
end

for d = 1:length(unique_wav) %detection loop
    file = unique_wav(d);
    dt_start = readDateTime(char((unique_wav(d)))); %start time of file, read in from filename
    PATH2WAV = fullfile(PATH2DATA,char((unique_wav(d))));
    [x] = audioread(PATH2WAV);%read in wav file
    [X] = audioread(PATH2WAV, 'native');
    p = utilities.amar2dB(X,5,-164);
    
    [M,q] = size(x); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds

    if plot_switch1 == 1
        figure(1)
        subplot(3,1,1)
        plot(t,x) %plot entire audio wav
        title("Entire audio wav")
        subplot(3,1,2)
        plot(t,X)
        title("Entire native audio wav")
        subplot(3,1,3)
        plot(t,p)
        title("SPL linear")
    end
    
   LowerStopbandFrequency = freq - bandpass_width-1000;
   LowerPassbandFrequency = freq - bandpass_width;
   UpperPassbandFrequency = freq + bandpass_width;
   UpperStopbandFrequency = freq + bandpass_width+1000;
        
   %%% create bandpass filter object if it doesn't exist already
   if isempty(bandpass_filter) || freq ~= filterFreq
      filterFreq = freq;
      bandpass_filter = designfilt(...
                 'bandpassfir',...
                 'StopbandFrequency1', LowerStopbandFrequency,...
                 'PassbandFrequency1', LowerPassbandFrequency,...
                 'PassbandFrequency2', UpperPassbandFrequency,...
                 'StopbandFrequency2', UpperStopbandFrequency,...
                 'StopbandAttenuation1', 60,...
                 'StopbandAttenuation2', 60,...
                 'PassbandRipple', 1,...
                 'DesignMethod', 'kaiserwin',...
                 'SampleRate', Fs...
                 );
   end 
   %Bandpass data
   x_freq = bandpass.noDelayFilt(bandpass_filter, x);
   
   Ms = max(abs(x_freq)); %gets maximum from bandpassed audio
   if Ms >= .95
       disp("WARNING: Audio saturated")
   end
        
   if plot_switch1 == 1
      figure(2)
      subplot(2,1,1)
      plot(x_freq) %plot bandpassed audio wav with mean removed
      title("Bandpassed audio wav")
   end
        
   if plot_switch1 == 1
      figure(2)
      subplot(2,1,2)
      spectrogram(x_freq,256,128,[],Fs,'yaxis')
      colorbar off
      ylim([LowerStopbandFrequency/1000 UpperStopbandFrequency/1000])
   end
   
   Filtered_peaks_wav = Filtered_PEAKS(Filtered_PEAKS.WavFiles == file,:);
   
   for ping = 1:height(Filtered_peaks_wav)
       ping_loc = Filtered_peaks_wav.peak_loc_freq(ping);
       ping_samps = Filtered_peaks_wav.FreqSAMPS90(ping);
       buffer = 512;
       ping_start_stop = [ping_loc-ceil(ping_samps/2), ping_loc + ceil(ping_samps/2)];
       ping_loc_clip = ceil((ping_samps/2)+1+buffer);
       ping_clip = [ping_start_stop(1)-buffer, ping_start_stop(2)+buffer];
       ping_x_clip = x_freq(ping_clip(1):ping_clip(2));
       
        if plot_switch1 == 1
            figure(2)
            subplot(2,1,1)
            hold on
            plot(ping_loc,max(x_freq),'*r')
            hold off
            figure(3)
            subplot(2,1,1)
            plot(ping_x_clip) %plot bandpassed audio wav with mean removed
            title("Bandpassed ping")
            hold on
            plot(ping_loc_clip,max(ping_x_clip),'*r')
            hold off
        end
        
        if plot_switch1 == 1
            figure(3)
            subplot(2,1,2)
            spectrogram(ping_x_clip,64,32,[],Fs,'yaxis')
            colorbar off
            ylim([LowerStopbandFrequency/1000 UpperStopbandFrequency/1000])
        end
        disp(' ')
        disp(file)
        ui = input('Enter 1 to validated ping, any other key to continue:','s');
        if ui == '1'
            Filtered_PEAKS.validated(ping) = 1;
            Filtered_PEAKS.reviewed(ping) = 1;
        elseif ui == 's'
        save(fullfile(PATH2OUTPUT,output_name), "Filtered_PEAKS")
        else
            Filtered_PEAKS.reviewed(ping) = 1;
            continue
        end
        
   end  % end ping loop
   
end                             %end filelist loop

save(fullfile(PATH2OUTPUT,output_name), "Filtered_PEAKS")

