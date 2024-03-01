%Extract_signal.m
%
% First attempt at using LTSA logger times to extract EK-60 pings to make
% standardized echo sounder ping

clear
close all

Fs = 256000;
%N=1024;                 % sample length
N=Fs;   %N = Fs pull of 1 second of data
N2=N/2;                 % step size = 30% overlap 30 seconds
frequencies =[18000, 38000, 50000, 70000, 120000];
bandpass_width = 1000; % +- width of bandpass filter

%threshold=dB2wav(110); % set threshold %%% need to create dB2wav function
%for AMAR



% load logs
PATH2LOG = 'D:\BW_ECHO_EXPERIMENT\COC_2020_09\COC_2020_09_echo_log.csv';
log = readtable(PATH2LOG);

%REMOVE FOR FUTURE USE
%Cludge to fix data directory restructure
log.InputFile = strrep(log.InputFile,"_ECHOSOUNDER_","_ECHO_");
%END CLUDGE

%load audio file
PATHfileList = log.InputFile;

%load strong ping
%P = audioread('D:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_STRONG_PING_TEMPLATE_1.wav');    %this in a main beam ping
P = audioread('D:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_REFLECTED_PING_TEMPLATE_1.wav'); %this in a reflected ping
MP = max(P);
%normalSP = SP/MSP(1); %normalize to max not nessecary if bandpassing
%according to Jinshan
%env_SP = envelope(normalSP,256000,'analytic'); %hilbert transform to extract envelope
%figure(8)
%plot(normalSP)
%ylim([-1.5 1.5])
%hold on
%plot(env_SP)
%hold off
b = conj(P(end:-1:1)); %inverse conjugate of the normalized ping
        
uniqueFileList = unique(PATHfileList); %get filelist from preliminary LTSA analysis

ping_indices =[];

for f = 1:length(uniqueFileList)%start filelist loop
%for f = 2 %this file has strong pings... some saturated though
    PATH2WAV = char(uniqueFileList(f));
    dt_start = readDateTime(PATH2WAV); %start time of file, read in from filename
    [x] = audioread(PATH2WAV); %read in wav file 
    disp(PATH2WAV);
    [M,q] = size(x);
    dt = 1/Fs;      %time between samples in seconds
    t = dt*(0:M-1)';%get time index in seconds
    x = detrend(x); %remove mean from audio %%%
    y1 = filter(b,1,x); %use match filter to filter audio for inverse conjugate of signal (echosounder ping)
    
    plot_switch1 = 0; %turns test plots on (1) or off (0)
    if plot_switch1 == 1
        figure(1)
        subplot(2,1,1);
        plot(x) %plot entire audio wav with mean removed and normalized to max
        title("Entire audio wav with mean removed and normalized to max")
        subplot(2,1,2);
        plot(y1) %plot match filtered audio
        title("Match filtered audio")
    end
    
    for freq = 1:length(frequencies) %Start frequency loop
    %for freq = 3
        freq_bins = [frequencies(freq)-bandpass_width frequencies(freq)+bandpass_width]; %create frequency bands using frequencies and width of bandpass
        x_freq = bandpass(x,freq_bins,Fs); %bandpass to isolate each frequency
        Ms = max(x_freq); %gets maximum from bandpassed audio
        if Ms == 1
            disp("WARNING: Audio saturated")
            
        end
        normalx = x_freq/Ms(1); %normalize bandpassed audio %needed?
        P_freq = bandpass(P,freq_bins,Fs); %bandpass signal (echosounder ping)
        MP_freq = max(P_freq); %maximum from template echosounder ping
        normalP_freq = P_freq/MP_freq; %normalize signal
        
        %%% attempt at moving search window loop %%%
         
        steps = [1:(length(normalx)/N2)]'; %define number of steps using number of samples in audio and user defined width of search window
        
        b_freq = conj(normalP_freq(end:-1:1)); %calculate inverse conjugate of bandpassed signal
        
        for n = 1:length(steps)-1
            nstart=(n-1)*N2+1;
            nend=nstart+N-1;
            window = nstart:nend;
            win_x = x(window); %select out window
            max_win = max(win_x);
            norm_win_x = win_x/max_win; %normalize window
            %figure(5)
            %subplot(2,1,1)
            %plot(norm_win_x)
            y2 = filter(b_freq,1,norm_win_x); %match filter
            %subplot(2,1,2)
            %plot(y2)
            %hold on
            fl1 = 60;
            [up1,lo1] = envelope(y2,fl1,'analytic');
            [Peak_val, peak_loc] = findpeaks(up1,'MinPeakProminence',8);
            %plot(up1)
            %hold on 
            if ~isempty(peak_loc)
               %plot(peak_loc,5,'r^') 
               [~,Max_ping] = max(Peak_val);
               max_ping_index = nstart+peak_loc(Max_ping);
               ping_indices = [ping_indices; max_ping_index];
               %hold off
               %pause() 
            else
               %hold off
               continue
            end
            %return
        end                     % end search window loop
        
   return     
        
        %%%%
        %env_normalSP_freq=envelope(normalSP,256000,'analytic'); %hilbert transform to extract envelope
        %figure(8)
        %plot(normalSP_freq)
        %ylim([-1.5 1.5])
        %hold on
        %plot(env_normalSP_freq)
        %hold off
        %%%
        %b_freq = conj(normalSP(end:-1:1)); %calculate the conjugate of the inverted signal
        %figure(2)
        %plot(normalx) %plot normalized wav with mean removed
        %hold on
        %y2 = filter(b_freq,1,x);
        %y = filter(b_freq,1,normalx); %match filter
       
        [Peak_val, peak_loc] = findpeaks(envelope(normalx),'MinPeakHeight',0.3,'MinPeakDistance',min_IPI);
        plot(peak_loc,0.5,'b^')
        hold on
        
        k = find(log.InputFile == uniqueFileList(f));
        Echo_times = log.StartTime(k,:);
        ping = seconds(Echo_times - dt_start);
        ping_samp = ping*Fs;
        plot(ping_samp,0,'r*')
        hold off
        for p = 1:length(ping) %Start ping loop
            %tping = ping(p);
            sample_ping = ping_samp(p);
            figure(1)
            hold on
            plot(sample_ping,0,"^")
            hold off
            sample_clip = [sample_ping-128000 sample_ping+128000];
            %clip = [tping-0.5 tping+0.5];
            xclip = x(ceil(min(sample_clip)):ceil(max(sample_clip))-1);
            test2 = findpeaks(xclip);
            %tclip = t(j);
            figure(2)
            plot(x);
            xlim(sample_clip)
            title(p + "/" + length(ping));
            hold on
            plot(128000,0,'r*')
            hold off
            [US,AAmp] = ginput;
            if ~isempty(US)
                for e = 1:length(US)  %Start user selected ping loop
                    sping = US(e);
                    clip = [sping-ser_win sping+ser_win];
                    sping_clip = x(ceil(min(clip)):ceil(max(clip))-1);
                    [test4, peaks_loc4] = findpeaks(envelope(sping_clip));
                    figure(3);
                    plot(x)
                    xlim(clip);
                    sping_env = envelope(sping_clip);
                    [m,i] = max(sping_env);
                    UserDiff = i-ser_win;
                    pingpeak = sping+UserDiff;
                    peak_clip = [pingpeak-ser_win pingpeak+ser_win];
                    ping_peak_x = x(ceil(min(peak_clip)):ceil(max(peak_clip))-1);
            
                    figure(4)
                    plot(ping_peak_x)
                    xlim([1 2*ser_win])
                    pause()
                end                   %End user selected ping loop
            end
         end                    %End ping loop
                                %clean up plots for next loop
        figure(1); clf;
        figure(2); clf;
        return
   end                         %end freq loop
end                             %end filelist loop

% click length ~0.002s


