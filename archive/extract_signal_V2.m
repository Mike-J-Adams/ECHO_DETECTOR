%Extract_signal.m
%
% First attempt at using LTSA logger times and matched filtering to extract EK-60 pings to make
% standardized echo sounder ping detector
tic
clear
close all

Fs = 256000;

%N=1024;                 % sample length
N=Fs;   %N = Fs pull of 1 seconds of data
N2=N/2;                 % step size = 50% overlap .5 seconds
frequencies =[18000, 38000, 50000, 70000, 120000];
bandpass_width = 2500; % +- width of bandpass filter

SNR_threshold = 10;
fl1= 60;

%set up Window
window = 0.004; %seconds
window_samps = window*Fs;

% load logs
PATH2LOG = 'E:\BW_ECHO_EXPERIMENT\COC_2020_09\COC_2020_09_echo_log.csv';
log = readtable(PATH2LOG);

%REMOVE FOR FUTURE USE
%Cludge to fix data directory restructure
log.InputFile = strrep(log.InputFile,"D:","E:");
log.InputFile = strrep(log.InputFile,"_ECHOSOUNDER_","_ECHO_");
%END CLUDGE

%load audio file
%PATHfileList = log.InputFile;
PATHfileList = dir("E:\BW_ECHO_EXPERIMENT\COC_2020_09\3DaySubset\*.wav");

%load strong ping
P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\Mean_Ping\MEAN_PING.wav');
%P = audioread('D:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_STRONG_PING_TEMPLATE_1.wav');    %this in a main beam ping
%P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_REFLECTED_PING_TEMPLATE_1.wav'); %this in a reflected ping
[MP, MPI] = max(P);
PT_window = MPI - ceil(window_samps/2):MPI + ceil(window_samps/2)-1;
PW = P(PT_window);
b = conj(PW(end:-1:1)); %inverse conjugate of the normalized ping all frequencies


%uniqueFileList = unique(PATHfileList);%get filelist from preliminary LTSA
%analysis
for f = 1:length(PATHfileList)%start filelist loop

%for f = 116 %this file has strong pings... some saturated though
    ping_indices =[];
    
    SNR_values = [];
    
    PINGS = [];
    SNR1 = [];
    SNR2 = [];
    RMS90 = [];
    DUR90 = [];
    SAMPS90 = [];
    
    %PATH2WAV = char(uniqueFileList(f));
    PATH2WAV = [PATHfileList(f).folder, '\', PATHfileList(f).name];
    dt_start = readDateTime(PATH2WAV); %start time of file, read in from filename
    [x] = audioread(PATH2WAV); %read in wav file 
    disp(PATH2WAV);
    [M,q] = size(x); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    %x = highpass(x,100,Fs); %takes too long...
    x = detrend(x); %remove mean from audio %%% takes a long time 21 seconds per file
    y1 = filter(b,1,x); %use match filter to filter audio for inverse conjugate of signal (echosounder ping)
    y1_abs = abs(y1);
    %[upall,loall] = envelope(y1_abs,fl1,'analytic');
    [Peak_val_all, peak_loc_all] = findpeaks(y1_abs,'MinPeakDistance',N);
   %[Peak_val_all, peak_loc_all] = findpeaks(upall);
    
    plot_switch1 = 0; %turns test plots on (1) or off (0)
    if plot_switch1 == 1
        figure(1)
        subplot(2,1,1);
        plot(t,x) %plot entire audio wav with mean removed 
        title("Entire audio wav with mean removed")
        subplot(2,1,2);
        plot(t,y1_abs) %plot match filtered audio
        hold on
        plot(peak_loc_all/Fs,max(y1_abs)*1.05,'^r')
        title("Match filtered audio")
        hold off
    end
 
    for n = 1:length(peak_loc_all)
    %for n = 1
            d_ping = peak_loc_all(n);
            ping_window = d_ping - ceil(window_samps/2):d_ping + ceil(window_samps/2)-1;
            noise_window = ping_window-(window_samps/2+1024);           
               if ping_window(end) > length(x)
                  continue
               end 
               if ping_window(1) < 1
                  disp("ping overlaps start of file")
                  continue
               end
               if noise_window(end) > length(x)
                   disp("noise overlaps end of file")
                   continue
               end 
               if noise_window(1) < 1
                  disp("noise overlaps start of file")
                  continue
               end
            ping = y1(ping_window);
            %%%%% extract 90% power
            cumPing = cumsum(abs(ping));
            maxPing = max(cumPing); %117
            minPing = min(cumPing);
            Ping951 = interp1(cumPing, cumPing, maxPing*0.95,'previous'); 
            Ping952 = interp1(cumPing, cumPing, maxPing*0.95,'nearest');  
            if Ping951 <= maxPing*0.95
                Ping95 = Ping952;
            elseif Ping951 > maxPing*0.95
                Ping95 = Ping951;
            end  
            Ping051 = interp1(cumPing, cumPing, maxPing*0.05,'next'); 
            Ping052 = interp1(cumPing, cumPing, maxPing*0.05,'nearest');  
            if Ping051 >= maxPing*0.05
                Ping05 = Ping052;
            elseif Ping051 < maxPing*0.05
                Ping05 = Ping052;
            end
            index95 = find(cumPing == Ping95);
            index05 = find(cumPing == Ping05);
            Ping90_window = [index05+min(ping_window):index95+min(ping_window)];
            Ping90 = y1(Ping90_window);
            rms90 = rms(Ping90);
            RMS90 = [RMS90, rms90];
            dur90 = length(Ping90)/Fs;
            DUR90 = [DUR90, dur90];
            samps90 = length(Ping90);
            SAMPS90 = [SAMPS90, samps90];
            env_ping = envelope(ping,10,'peak');
            max_ping = max(abs(ping));
            noise_sample = y1(noise_window);
            SNR = rssq(ping(:))/rssq(noise_sample(:));
            SNR_values = [SNR_values; SNR];
            noise90_window = Ping90_window - (window_samps/2+1024);
            noise90 = y1(noise90_window);
            snr1 = snr(Ping90,noise90);
            SNR1 = [SNR1, snr1];
            snr2 = snr(ping,noise_sample);
            SNR2 = [SNR2, snr2];
            if plot_switch1 == 1
                figure(2)
                subplot(2,1,1);
                plot(Ping90_window,Ping90) %plot matched filtered audio 
                title("Ping")
                subplot(2,1,2);
                plot(noise90_window,noise90) %plot match filtered audio
                title("Noise sample")
            end
            
            if SNR >= SNR_threshold
               %[~,Max_ping] = max(Peak_val); %takes max ping
               ping_index = d_ping; %finds index in data
               ping_indices = [ping_indices; ping_index];
               
            end
           
    end
    ping_indices = unique(ping_indices);
    ping_clips = [ping_indices - ceil(window_samps/2) ping_indices + ceil(window_samps/2)];
    
    
   
    for freq = 1:length(frequencies) %Start frequency loop
    %for freq = 1
        freq_bins = [frequencies(freq)-bandpass_width frequencies(freq)+bandpass_width]; %create frequency bands using frequencies and width of bandpass
        %%%% Think about bandpass more...
        x_freq = bandpass(x,freq_bins,Fs); %bandpass to isolate each frequency
        %%%%
        
        Ms = max(x_freq); %gets maximum from bandpassed audio
        if Ms >= .95
            disp("WARNING: Audio saturated")
        end
        normalx = x_freq/Ms(1); %normalize bandpassed audio %needed?
        P_freq = bandpass(PW,freq_bins,Fs); %bandpass signal (echosounder ping)
        MP_freq = max(P_freq); %maximum from template echosounder ping
        normalP_freq = P_freq/MP_freq; %normalize signal

    %%% attempt at moving search window loop %%%
        steps = [1:ceil((length(normalx)/N2))]'; %define number of steps using number of samples in audio and user defined width of search window
    %%%     
        b_freq = conj(normalP_freq(end:-1:1)); %calculate inverse conjugate of bandpassed signal
        fl1 = 60;
        y2 = filter(b_freq,1,normalx); %match filter

        if plot_switch1 == 1
            figure(3)
            subplot(2,1,1);
            plot(normalx) %plot entire audio wav with mean removed and normalized to max
            title("Bandpassed audio wav with mean removed and normalized")
            subplot(2,1,2);
            plot(y2) %plot match filtered audio
            title("Match filtered bandpassed audio")
        end
        
        y2_abs = abs(y2);
        
        %[upfreq,lofreq] = envelope(y2_abs,fl1,'analytic');
        [Peak_val_freq, peak_loc_freq] = findpeaks(y2_abs,'MinPeakDistance',N);
        %[Peak_val_all, peak_loc_all] = findpeaks(upall);
        
        SNR_values_freq = [];
        FreqPINGS = [];
        FreqSNR1 = [];
        FreqSNR2 = [];
        FreqRMS90 = [];
        FreqDUR90 = [];
        FreqSAMPS90 = [];
        ping_indices_freq = [];
        
        for n = 1:length(peak_loc_freq)
            d_ping_freq = peak_loc_freq(n);
            ping_window_freq = d_ping_freq - ceil(window_samps/2):d_ping_freq + ceil(window_samps/2)-1;
            noise_window_freq = ping_window_freq-(window_samps/2+1024);           
               if ping_window_freq(end) > length(y2_abs)
                  continue
               end 
               if ping_window_freq(1) < 1
                  continue
               end
               if noise_window_freq(end) > length(y2_abs)
                   continue
               end 
               if noise_window_freq(1) < 1
                  continue
               end
            ping_freq = y2_abs(ping_window_freq);
            max_ping_freq = max(abs(ping_freq));
            noise_sample_freq = y2_abs(noise_window_freq);
            SNR_freq = rssq(ping_freq(:))/rssq(noise_sample_freq(:));
            SNR_values_freq = [SNR_values_freq; SNR_freq];
            
            if plot_switch1 == 1
                figure(4)
                subplot(2,1,1);
                plot(ping_freq) %plot matched filtered audio 
                title("Ping")
                subplot(2,1,2);
                plot(noise_sample_freq) %plot match filtered audio
                title("Noise sample")
            end
            
            cum_FreqPing = cumsum(abs(ping_freq));
            maxFreqPing = max(cum_FreqPing);
            minFreqPing = min(cum_FreqPing);
            
            %%% This part is confusing
            FreqPing951 = interp1(cum_FreqPing, cum_FreqPing, maxFreqPing*0.95,'previous'); 
            FreqPing952 = interp1(cum_FreqPing, cum_FreqPing, maxFreqPing*0.95,'nearest');  
            if FreqPing951 <= maxFreqPing*0.95
                FreqPing95 = FreqPing952;
            elseif FreqPing951 > maxFreqPing*0.95
                FreqPing95 = FreqPing951;
            end  
            FreqPing051 = interp1(cum_FreqPing, cum_FreqPing, maxFreqPing*0.05,'next'); 
            FreqPing052 = interp1(cum_FreqPing, cum_FreqPing, maxFreqPing*0.05,'nearest');  
            if FreqPing051 <= maxFreqPing*0.05
                FreqPing05 = FreqPing052;
            elseif FreqPing051 > maxFreqPing*0.05
                FreqPing05 = FreqPing051;
            end
            %%%
            
            Freqindex95 = find(cum_FreqPing == FreqPing95);
            Freqindex05 = find(cum_FreqPing == FreqPing05);
            FreqPing90_window = [Freqindex05+min(ping_window_freq):Freqindex95+min(ping_window_freq)];
            FreqPing90 = y2(FreqPing90_window);
            Freqrms90 = rms(FreqPing90);
            FreqRMS90 = [FreqRMS90, Freqrms90];
            Freqdur90 = length(FreqPing90)/Fs;
            FreqDUR90 = [FreqDUR90, Freqdur90];
            Freqsamps90 = length(FreqPing90);
            FreqSAMPS90 = [FreqSAMPS90, Freqsamps90];
            Freqenv_ping = envelope(ping_freq,10,'peak');
            Freqmax_ping = max(abs(ping_freq));
            %Freq_noise_sample = y1(Freqnoise_window);
            %FreqSNR = rssq(ping_freq(:))/rssq(Freq_noise_sample(:));
            %SNR_values = [SNR_values; SNR];
            Freqnoise90_window = FreqPing90_window - (window_samps/2+1024);
            Freqnoise90 = y2(Freqnoise90_window);
            Freqsnr1 = snr(FreqPing90,Freqnoise90);
            FreqSNR1 = [FreqSNR1, Freqsnr1];
            Freqsnr2 = snr(ping_freq,noise_sample_freq);
            FreqSNR2 = [FreqSNR2, Freqsnr2];
            
            if plot_switch1 == 1
                figure(5)
                subplot(2,1,1);
                plot(FreqPing90_window,FreqPing90) %plot matched filtered audio 
                title("Ping90")
                subplot(2,1,2);
                plot(Freqnoise90_window,Freqnoise90) %plot match filtered audio
                title("Noise sample90")
            end
            
            
            if Freqsnr1 >= SNR_threshold
               %[~,Max_ping] = max(Peak_val); %takes max ping
               ping_index_freq = d_ping_freq; %finds index in data
               ping_indices_freq = [ping_indices_freq; ping_index_freq];
            else
               disp("No ping")
            end
         end
    ping_indices_freq = unique(ping_indices_freq);
    ping_clips_freq = [ping_indices_freq - ceil(window_samps/2) ping_indices_freq + ceil(window_samps/2)];
             
    freq_detect = table(peak_loc_freq,FreqDUR90',FreqRMS90',FreqSAMPS90',FreqSNR1');
    
    freq_detect_SNR = freq_detect(freq_detect.Var5 >=10,:);
           % figure(4)
            %subplot(2,1,1)
            %plot(ping)
            %subplot(2,1,2)
            %plot(noise_sample)
      
                                        %clean up plots for next loop

            if ~isempty(freq_detect_SNR)
                %figure(1)
                %subplot(2,1,2)
                %hold on
                %plot(ping_indices/Fs,max(y1_abs)*1.10,'^g')
                %hold off
                slashIdx = strfind(PATH2WAV, '\'); 
                pathSegment = PATH2WAV(slashIdx(end)+1:end);
                file_n = split(pathSegment,'.');
                new_freq_file = ['match_' num2str(frequencies(freq)) '_' char(file_n(1)) '.' char(file_n(2)) '.mat'];
                %new_file = ['match_ALL_'  char(file_n(1)) '.' char(file_n(2)) '.mat'];
                new_table = ['match_freq_' num2str(frequencies(freq)) '_' char(file_n(1)) '.' char(file_n(2)) '.mat'];
                Path2Output = [PATH2WAV(1:slashIdx(end-1)) 'OUTPUT\'];
                if ~exist(Path2Output, 'dir')  
                    mkdir(Path2Output);
                end
                save_path1 = string([Path2Output new_freq_file]);
                save_path2 = string([Path2Output new_table]);
                save(save_path1, 'ping_indices')
                save(save_path2, 'freq_detect_SNR')
            end
            
            figure(1); clf;
            figure(2); clf;
   end%end freq loop
end   %end filelist loop

% click length ~0.002s

toc
