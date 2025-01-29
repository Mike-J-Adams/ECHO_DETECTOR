%fill_gaps.m

%fill gaps in detections

clear
close all

import('utilities.readDateTime');

Fs = 256000;
%N=1024;                % sample length
N=Fs*30;                 %N = Fs pull of 10 second of data
N2=N/2;                 
frequencies =[18000, 38000, 50000, 70000, 120000];
bandpass_width = 5000;  % +- width of bandpass filter

%%% Change to review single frequency, leave empty to review all
freq = 18000;
%%%

PATH2OUTPUT = 'F:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\COC_2020_09';
PATH2DETECTIONS = 'F:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\COC_2020_09\GBK_EK60_DETECTIONS_FILTERED_120kHz_VALIDATED.mat';
PATH2DATA = 'F:\BW_ECHO_EXPERIMENT\WAV_FILES\COC_2020_09\3DAYSUBSET';
PATH2ADDEDPINGS = 'F:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\COC_2020_09\ADDED_18kHz_PINGS.mat';
load(PATH2DETECTIONS);
%filter out already validated detections
Filtered_peaks_wav_reviewed = Filtered_peaks_wav_reviewed(Filtered_peaks_wav_reviewed.validated == 1,:);

if isfile(PATH2ADDEDPINGS)
   load(PATH2ADDEDPINGS); 
end

%%% create empty variable to store bandpass filter object
bandpass_filter = [];
filterFreq = [];

unique_wav = unique(Filtered_peaks_wav_reviewed.WavFiles);

new_pings = [];

for d = 5:length(unique_wav) %detection loop
    file = unique_wav(d);
    dt_start = readDateTime(char((unique_wav(d)))); %start time of file, read in from filename
    PATH2WAV = fullfile(PATH2DATA,char((unique_wav(d))));
    [x] = audioread(PATH2WAV);%read in wav file
    [X] = audioread(PATH2WAV, 'native');
    ainfo = audioinfo(PATH2WAV);
    X_adjusted = X/256;
    
    disp(' ')
    disp(file)
      
    [M,q] = size(x); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds

   figure(1)
   subplot(2,1,1)
   plot(t,x) %plot entire audio wav
   title("Entire audio wav")
   subplot(2,1,2)
   plot(t,X)
   title("Entire native audio wav")
       
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

   Filtered_peaks_wav = Filtered_peaks_wav_reviewed(Filtered_peaks_wav_reviewed.WavFiles == file,:);
   
    
   figure(2)
   subplot(2,1,1)
   plot(x_freq) %plot bandpassed audio wav with mean removed
   title("Bandpassed audio wav")
   hold on 
   plot(Filtered_peaks_wav.adjusted_ping_loc,max(x_freq),'*r')
   hold off
   subplot(2,1,2)
   spectrogram(x_freq,256,128,[],Fs,'yaxis')
   colorbar off
   ylim([LowerStopbandFrequency/1000 UpperStopbandFrequency/1000])


   
   figure(3)
   %subplot(2,1,1)
   plot(x_freq) %plot bandpassed audio wav with mean removed
   title("Bandpassed ping wav")
   hold on 
   plot(Filtered_peaks_wav.adjusted_ping_loc,max(x_freq),'*r')
   hold off
   
   figure(3)
   [filli,fillv] = ginput;
   
   figure(4)
   plot(x_freq)
   
   figure(5)
   plot(x_freq)
   
   
   %for ping = 16:height(Filtered_peaks_wav)
   for ping = 1:length(filli)
       ping_loc = filli(ping);
       ping_start_stop = [ping_loc-N2, ping_loc+N2];
       ping_loc_clip = [];
       
       if ping_start_stop(1) <= 0
        ping_start_stop(1) = 1;
       end
       if ping_start_stop(2) > ainfo.TotalSamples
        ping_start_stop(2) = ainfo.TotalSamples;
       end
       
       ping_x_clip = x_freq(ping_start_stop(1):ping_start_stop(2)-1);
       time_x_clip = t(ping_start_stop(1):ping_start_stop(2)-1);
       
       %figure(2)
       %subplot(2,1,1)
       %hold on 
       %plot(ping_loc,max(x_freq),'*g')
       %hold off
       
       figure(3)
       %subplot(2,1,1)
       hold on 
       %plot(ping_loc,max(x_freq),'*g')
       xlim([ping_start_stop(1) ping_start_stop(2)])
       hold off
       %figure(3)
       %subplot(2,1,1)
       %plot(time_x_clip, ping_x_clip) %plot bandpassed audio wav with mean removed
       %plot(x_freq) %plot bandpassed audio wav with mean removed
       %title("Bandpassed ping")
       
              
       %figure(3)
       %subplot(2,1,2)
       %spectrogram(ping_x_clip,8192,8192/2,[],Fs,'yaxis')
       %colorbar off
       %ylim([LowerStopbandFrequency/1000 UpperStopbandFrequency/1000])
       

       
       figure(3)
       %subplot(2,1,1)
       [pin,pvn] = ginput;
       window =Fs*0.25;
       
       if isempty(pin)
           continue
       end
       
       
       
       
       for pindex = 1:length(pin)
          subclip_start_stop = [floor(pin(pindex) - window/2), ceil(pin(pindex)+window/2)-1];
          
          if subclip_start_stop(1) <= 0
            subclip_start_stop(1) = 1;
          end
          if subclip_start_stop(2) > length(x_freq)
            subclip_start_stop(2) = length(x_freq);
          end
          subclip = x_freq(subclip_start_stop(1):subclip_start_stop(2)-1);
          [Msub,qsub] = size(subclip); %get size length of audio
          tsub = dt*(0:Msub-1)';%get time index in seconds
          
          figure(4)
          xlim([subclip_start_stop(1),subclip_start_stop(2)])
          
          [max_added_V,max_added_I] = max(x_freq(subclip_start_stop(1):subclip_start_stop(2)));
          new_ping_loc_plot = subclip_start_stop(1) + max_added_I;
          new_ping_loc(pindex) = subclip_start_stop(1) + max_added_I;
          
          hold on
          plot(new_ping_loc_plot,max_added_V,'*r')
          hold off
          
          addedping_start_stop = [new_ping_loc_plot - window/16, new_ping_loc_plot + window/16];
          
          if addedping_start_stop(1) <= 0
            addedping_start_stop(1) = 1;
          end
          if addedping_start_stop(2) > length(x_freq)
            addedping_start_stop(2) = length(x_freq);
          end
          addedping = x_freq(addedping_start_stop(1):addedping_start_stop(2)-1);
          [Mnew,qnew] = size(addedping); %get size length of audio
          tnew = dt*(0:Mnew-1)';%get time index in seconds
                 
          figure(5)
          xlim([addedping_start_stop(1),addedping_start_stop(2)])
          hold on
          plot(new_ping_loc_plot,max_added_V,'*r')
          hold off
          
          ui = input('Enter 1 to validated ping, any other key to continue:','s');
          if ui == '1'
              validated = 1;
              figure(3)
              hold on
              plot(new_ping_loc_plot,max(x_freq),'*g')
              hold off
          else
              validated = 0;
          end
        temp_new_pings = [file new_ping_loc_plot validated freq];
        new_pings = [new_pings;temp_new_pings];
       end
       new_pings = new_pings(str2double(new_pings(:,3))==1,:);
       output_file = char(fullfile(PATH2OUTPUT,['ADDED_',num2str(freq/1000),'kHz_PINGS']));
       save(output_file,'new_pings')
       
   end  
   

   
end                             %end filelist loop



