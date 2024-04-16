%Extract_signal.m
%
% First attempt at using matched filtering to find EK-60 pings in acoustic
% data

tic
clear
close all


%set channel
channel = 1;
%set sample rate

Fs = 256000;

N=Fs*3;   %N = Fs pull 1 second of data %%% currently set at 3 seconds

% set frequencies of interest
frequencies =[18000, 38000, 50000, 70000, 120000]; %ek60
% +- width of bandpass filter
bandpass_width = 2500; 

%set up sample extraction window
window = 0.003; %seconds
window_samps = window*Fs;

%Set path to .wav files
PATH2DATA = "E:\BW_ECHO_EXPERIMENT\COC_2020_09\3DaySubset\"; %ek60 all files
%PATH2Data = "E:\BW_ECHO_EXPERIMENT\COC_2020_09\strong_pings_subset\*.wav"; %bigelow around
PATHfileList = dir(fullfile(PATH2DATA,"*.wav"));
%Set path to output
PATH2OUTPUT = "E:\BW_ECHO_EXPERIMENT\COC_2020_09\OUTPUT2"; 

%turns debugging plots on (1) or off (0)
plot_switch1 = 0; 

%restart logic
[iStart] = utilities.restart_logic(PATH2OUTPUT,PATH2DATA);

%load template mean ping to be used in matched filter
P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\Mean_Ping\MEAN_PING.wav'); %ek60 ping
%P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\Mean_Ping\18kHzPing.wav'); %ek60 ping pure 18kHz sine
%P = audioread('D:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_STRONG_PING_TEMPLATE_1.wav');    %this in a main beam ping
%P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_REFLECTED_PING_TEMPLATE_1.wav'); %this in a reflected ping
[MP, MPI] = max(abs(P)); %get maximum amplitude

%Clean up template duration to match window size
PT_window = MPI - ceil(window_samps/2):MPI + ceil(window_samps/2)-1;
PW = P(PT_window);

%%% create empty variable to store bandpass filter object
bandpass_filter = [];
filterFreq = [];

%analysis
for f = iStart:length(PATHfileList)%start filelist loop
    
    PATH2WAV = [PATHfileList(f).folder, '\', PATHfileList(f).name];
    %start time of file, read in from filename
    dt_start = utilities.readDateTime(PATH2WAV);
    %read in wav file
    [x] = audioread(PATH2WAV);  
    disp(PATH2WAV);

    [M,q] = size(x); %get size length of audio
    if q > 1
        x = x(:,channel);
    end
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    %x = highpass(x,100,Fs); %takes too long...
    %x = detrend(x); %remove mean from audio %%% takes a long time 21
    %seconds per file. I think this was dumb...
    
    plot_switch1 = 0; %turns test plots on (1) or off (0)

    for freq = 1:length(frequencies) %Start frequency loop
    %for freq = 1
        %define stop and pass bands for bandpass filter
        LowerStopbandFrequency = frequencies(freq) - bandpass_width-1000;
        LowerPassbandFrequency = frequencies(freq) - bandpass_width;
        UpperPassbandFrequency = frequencies(freq) + bandpass_width;
        UpperStopbandFrequency = frequencies(freq) + bandpass_width+1000;
        
        %%% create bandpass filter object if it doesn't exist already
        if isempty(bandpass_filter) || frequencies(freq) ~= filterFreq
           filterFreq = frequencies(freq);
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
        x_BP = bandpass.noDelayFilt(bandpass_filter, x);
        %get maximum from bandpassed audio and check if audio is saturated
        Ms = max(x_BP); 
        if Ms >= .95
            disp("WARNING: Audio saturated")
        end
        
        %normalize bandpassed audio %needed?
        normalx = x_BP/Ms(1);
        
        %bandpass signal (echosounder ping)
        P_freq = bandpass.noDelayFilt(bandpass_filter,PW);
        %maximum from template echosounder ping
        MP_freq = max(P_freq); 
        %normalize signal
        normalP_freq = P_freq/MP_freq; 
        %calculate inverse conjugate of bandpassed signal
        b_freq = conj(normalP_freq(end:-1:1)); 
       
        %match filter
        y2 = filter(b_freq,1,normalx); 

        if plot_switch1 == 1
            figure(1)
            subplot(2,1,1);
            plot(normalx) %plot entire audio wav with mean removed and normalized to max
            title("Bandpassed audio wav with mean removed and normalized")
            subplot(2,1,2);
            plot(y2) %plot match filtered audio
            title("Match filtered bandpassed audio")
        end
        
        %approximate energy by squaring amplitude
        y2_2 = y2.^2;
        
        % Find peaks in squared matched filter response greater than 3 second apart 
        [Peak_val_freq, peak_loc_freq] = findpeaks(y2_2,'MinPeakDistance',N);
        %potentially split script here for ease of use...
 
        
        
        
        %initialize storage containers
        peaks = table(peak_loc_freq, Peak_val_freq);
        peaks.retained = zeros(size(peak_loc_freq));
        peaks.FreqSNR1 = zeros(size(peak_loc_freq));
        peaks.FreqSNR2 = zeros(size(peak_loc_freq));
        peaks.FreqSNR3 = zeros(size(peak_loc_freq));
        peaks.FreqRMS90 = zeros(size(peak_loc_freq));
        peaks.FreqDUR90 = zeros(size(peak_loc_freq));
        peaks.FreqSAMPS90 = zeros(size(peak_loc_freq));
        
        %start looping through peaks 
        for n = 1:length(peak_loc_freq)
            %get peak index
            d_ping_freq = peak_loc_freq(n);
            %define window to extract rough peak 
            ping_window_freq = d_ping_freq - ceil(window_samps/2):d_ping_freq + ceil(window_samps/2)-1;
           
            if ping_window_freq(end) > length(y2_2)
               continue
            end 
            if ping_window_freq(1) < 1
               continue
            end
            
            peaks.retained(n) = 1;
            ping_freq = y2(ping_window_freq);
            max_ping_freq = max(abs(ping_freq));
            
            if plot_switch1 == 1
                figure(4)
                plot(ping_freq) %plot matched filtered audio 
                title("Ping")
            end
            
            [Start90, Stop90] = utilities.calcEng(ping_freq,90);
            RelativeStart90 = Start90+min(ping_window_freq);
            RelativeStop90 = Stop90+min(ping_window_freq);
            
            if plot_switch1 == 1
               figure(4)
               hold on
               xline(Start90)
               xline(Stop90)
               hold off
            end            
            
            [FreqPing90, Freqnoise90] = snr.extractSN(x, Fs, RelativeStart90, RelativeStop90, 1024, 512, bandpass_filter, 'samples');
            [snr_dB, snr_adjusted_dB] = snr.calculateSNR(FreqPing90, Freqnoise90);
            if isempty(snr_dB)
                snr_dB = NaN;
                snr_adjusted_dB = NaN;
            end
            
            Freqrms90 = rms(FreqPing90);
            peaks.FreqRMS90(n) = Freqrms90;
            Freqdur90 = length(FreqPing90)/Fs;
            peaks.FreqDUR90(n) = Freqdur90;
            Freqsamps90 = length(FreqPing90);
            peaks.FreqSAMPS90(n) = Freqsamps90;
            Freqenv_ping = envelope(ping_freq,10,'peak');
            Freqmax_ping = max(abs(ping_freq));
            
            Freqsnr1 = snr(FreqPing90,Freqnoise90);
            peaks.FreqSNR1(n) = Freqsnr1;
            
            peaks.FreqSNR2(n) = snr_dB;
            peaks.FreqSNR3(n) = snr_adjusted_dB;
            
            if plot_switch1 == 1
                figure(5)
                subplot(2,1,1);
                plot(FreqPing90)  
                title("Ping90")
                subplot(2,1,2);
                plot(Freqnoise90) %plot match filtered audio
                title("Noise sample90")
            end
        end   %end peak loop         
              % Save detection output!                      
               slashIdx = strfind(PATH2WAV, '\'); 
               pathSegment = PATH2WAV(slashIdx(end)+1:end);
               file_n = split(pathSegment,'.');
               new_table = ['match_table_freq_' num2str(frequencies(freq)) '_' char(file_n(1)) '.' char(file_n(2)) '.mat'];
               Path2Output = [PATH2WAV(1:slashIdx(end-1)) 'OUTPUT2\'];
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



toc
