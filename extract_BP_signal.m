%Extract_signal.m
%
% First attempt at using matched filtering to find EK-60 pings in acoustic
% data

tic
clear
close all

Fs = 256000;

%N=1024;                 % sample length
N=Fs;   %N = Fs pull 1 second of data
N2=N/2;                 % step size = 50% overlap .5 seconds
frequencies =[18000, 38000, 50000, 70000, 120000]; %ek60
bandpass_width = 2500; % +- width of bandpass filter

SNR_threshold = 10;
fl1= 60;

%set up Window
window = 0.004; %seconds
window_samps = window*Fs;

%load audio file
%PATHfileList = log.InputFile;
PATH2DATA = "E:\BW_ECHO_EXPERIMENT\COC_2020_09\3DaySubset\"; %ek60 all files
%PATH2Data = "E:\BW_ECHO_EXPERIMENT\COC_2020_09\strong_pings_subset\*.wav"; %bigelow around
PATHfileList = dir(PATH2DATA);

%restart logic
PATH2OUTPUT = "E:\BW_ECHO_EXPERIMENT\COC_2020_09\OUTPUT2"; 

[iStart] = utilities.restart_logic(PATH2OUTPUT,PATH2DATA);

%load template mean ping
P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\Mean_Ping\MEAN_PING.wav'); %ek60 ping
%P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\Mean_Ping\18kHzPing.wav'); %ek60 ping pure 18kHz sine
%P = audioread('D:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_STRONG_PING_TEMPLATE_1.wav');    %this in a main beam ping
%P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_REFLECTED_PING_TEMPLATE_1.wav'); %this in a reflected ping
[MP, MPI] = max(abs(P)); %get maximum amplitude  
%Clean up template duration to match window size
PT_window = MPI - ceil(window_samps/2):MPI + ceil(window_samps/2)-1;
PW = P(PT_window);
b = conj(PW(end:-1:1)); %inverse conjugate of the normalized ping all frequencies

%analysis
for f = iStart:length(PATHfileList)%start filelist loop
    
    PATH2WAV = [PATHfileList(f).folder, '\', PATHfileList(f).name];
    dt_start = utilities.readDateTime(PATH2WAV); %start time of file, read in from filename
    [x] = audioread(PATH2WAV); %read in wav file 
    disp(PATH2WAV);
    [M,q] = size(x); %get size length of audio
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    %x = highpass(x,100,Fs); %takes too long...
    %x = detrend(x); %remove mean from audio %%% takes a long time 21
    %seconds per file. I think this was dumb...
    
    %%% create empty variable to store bandpass filter object
    bandpass_filter = [];
    
    plot_switch1 = 0; %turns test plots on (1) or off (0)
      
    for freq = 1:length(frequencies) %Start frequency loop
    %for freq = 1
        freq_bins = [frequencies(freq)-bandpass_width frequencies(freq)+bandpass_width]; %create frequency bands using frequencies and width of bandpass
        %%%% Think about bandpass more...
        x_BP = bandpass(x,freq_bins,Fs); %bandpass to isolate each frequency
        %%%%
        Ms = max(x_BP); %gets maximum from bandpassed audio
        if Ms >= .95
            disp("WARNING: Audio saturated")
        end
        normalx = x_BP/Ms(1); %normalize bandpassed audio %needed?
        P_freq = bandpass(PW,freq_bins,Fs); %bandpass signal (echosounder ping)
        MP_freq = max(P_freq); %maximum from template echosounder ping
        normalP_freq = P_freq/MP_freq; %normalize signal
    
        b_freq = conj(normalP_freq(end:-1:1)); %calculate inverse conjugate of bandpassed signal
        fl1 = 60;
        y2 = filter(b_freq,1,normalx); %match filter

        if plot_switch1 == 1
            figure(1)
            subplot(2,1,1);
            plot(normalx) %plot entire audio wav with mean removed and normalized to max
            title("Bandpassed audio wav with mean removed and normalized")
            subplot(2,1,2);
            plot(y2) %plot match filtered audio
            title("Match filtered bandpassed audio")
        end
        
        y2_2 = y2.^2; %approximate energy by squaring amplitude
        

        [Peak_val_freq, peak_loc_freq] = findpeaks(y2_2,'MinPeakDistance',N);

        %initialize storage containers
        peaks = table(peak_loc_freq, Peak_val_freq);
        peaks.retained = zeros(size(peak_loc_freq));
        peaks.FreqSNR1 = zeros(size(peak_loc_freq));
        peaks.FreqSNR2 = zeros(size(peak_loc_freq));
        peaks.FreqRMS90 = zeros(size(peak_loc_freq));
        peaks.FreqDUR90 = zeros(size(peak_loc_freq));
        peaks.FreqSAMPS90 = zeros(size(peak_loc_freq));

        for n = 1:length(peak_loc_freq)
            d_ping_freq = peak_loc_freq(n);
            ping_window_freq = d_ping_freq - ceil(window_samps/2):d_ping_freq + ceil(window_samps/2)-1;
            noise_window_freq = ping_window_freq-(window_samps/2+1024);           
               if ping_window_freq(end) > length(y2_2)
                  continue
               end 
               if ping_window_freq(1) < 1
                  continue
               end
               if noise_window_freq(end) > length(y2_2)
                   continue
               end 
               if noise_window_freq(1) < 1
                  continue
               end
            peaks.retained(n) = 1;
            ping_freq = y2_2(ping_window_freq);
            max_ping_freq = max(abs(ping_freq));
            noise_sample_freq = y2_2(noise_window_freq);
            
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
            
            %%% 
            max95Ping = maxFreqPing*0.95;
            max05Ping = maxFreqPing*0.05;
            
            %find closest value to 95%
            V = max95Ping;
            C = cum_FreqPing;
            A = repmat(C,[1 length(V)]);
            [~,closestIndex] = min(abs(A-V'));
            FreqPing95 = C(closestIndex); 
            
            %find closest value to 5%
            V = max05Ping;
            C = cum_FreqPing;
            A = repmat(C,[1 length(V)]);
            [~,closestIndex] = min(abs(A-V'));
            FreqPing05 = C(closestIndex); 
            %%%
            
            Freqindex95 = find(cum_FreqPing == FreqPing95);
            Freqindex05 = find(cum_FreqPing == FreqPing05);
            FreqPing90_window = [Freqindex05+min(ping_window_freq):Freqindex95+min(ping_window_freq)];
            Start90 = Freqindex05+min(ping_window_freq);
            Stop90 = Freqindex95+min(ping_window_freq);
            FreqPing90 = y2(FreqPing90_window);
                       
            Freqrms90 = rms(FreqPing90);
            peaks.FreqRMS90(n) = Freqrms90;
            Freqdur90 = length(FreqPing90)/Fs;
            peaks.FreqDUR90(n) = Freqdur90;
            Freqsamps90 = length(FreqPing90);
            peaks.FreqSAMPS90(n) = Freqsamps90;
            Freqenv_ping = envelope(ping_freq,10,'peak');
            Freqmax_ping = max(abs(ping_freq));
            Freqnoise90_window = FreqPing90_window - (window_samps/2+1024);
            Freqnoise90 = y2(Freqnoise90_window);
            Freqsnr1 = snr(FreqPing90,Freqnoise90);
            peaks.FreqSNR1(n) = Freqsnr1;
            Freqsnr2 = snr(ping_freq,noise_sample_freq);
            peaks.FreqSNR2(n) = Freqsnr2;
            
            if plot_switch1 == 1
                figure(5)
                subplot(2,1,1);
                plot(FreqPing90_window,FreqPing90) %plot matched filtered audio 
                title("Ping90")
                subplot(2,1,2);
                plot(Freqnoise90_window,Freqnoise90) %plot match filtered audio
                title("Noise sample90")
            end
        end   %end peak loop         
                                    
               slashIdx = strfind(PATH2WAV, '\'); 
               pathSegment = PATH2WAV(slashIdx(end)+1:end);
               file_n = split(pathSegment,'.');
               new_table = ['match_table_freq_' num2str(frequencies(freq)) '_' char(file_n(1)) '.' char(file_n(2)) '.mat'];
               Path2Output = [PATH2WAV(1:slashIdx(end-1)) 'OUTPUT3\'];
               if ~exist(Path2Output, 'dir')  
                  mkdir(Path2Output);
               end
               save_path2 = string([Path2Output new_table]);
               save(save_path2, 'peaks') 
               %clean up plots for next loo
            
               %figure(1); clf;
               %figure(2); clf;
   end%end freq loop
end   %end filelist loop

% click length ~0.002s

toc
