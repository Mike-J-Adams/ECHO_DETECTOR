%Extract_signal.m
%
% First attempt at using LTSA logger times to extract EK-60 pings to make
% standardized echo sounder ping

clear
close all

Fs = 256000;
%N=1024;                 % sample length
N=Fs;   %N = Fs pull of 1 second of data
N2=N/2;                 % step size = 50% overlap .5 seconds
frequencies =[18000, 38000, 50000, 70000, 120000];
bandpass_width = 1000; % +- width of bandpass filter

%threshold=dB2wav(110); % set threshold %%% need to create dB2wav function
%for AMAR



% load logs
PATH2LOG = 'E:\BW_ECHO_EXPERIMENT\COC_2020_09\COC_2020_09_echo_log.csv';
log = readtable(PATH2LOG);

%REMOVE FOR FUTURE USE
%Cludge to fix data directory restructure
log.InputFile = strrep(log.InputFile,"D:","E:");
log.InputFile = strrep(log.InputFile,"_ECHOSOUNDER_","_ECHO_");
%END CLUDGE

%load audio file
PATHfileList = log.InputFile;

%load strong ping
%P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_STRONG_PING_TEMPLATE_1.wav');    %this in a main beam ping
P = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\AMAR538_REFLECTED_PING_TEMPLATE_1.wav'); %this in a reflected ping
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


%ping_indices =[];

%for f = 1:length(uniqueFileList)%start filelist loop
for f = 2 %this file has strong pings... some saturated though
    ping_indices =[];
    PATH2WAV = char(uniqueFileList(f));
    dt_start = readDateTime(PATH2WAV); %start time of file, read in from filename
    [x] = audioread(PATH2WAV); %read in wav file 
    disp(PATH2WAV);
    [M,q] = size(x); %get size length of audio
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
    
    %for freq = 1:length(frequencies) %Start frequency loop
    for freq = 3
        freq_bins = [frequencies(freq)-bandpass_width frequencies(freq)+bandpass_width]; %create frequency bands using frequencies and width of bandpass
        x_freq = bandpass(x,freq_bins,Fs); %bandpass to isolate each frequency
        Ms = max(x_freq); %gets maximum from bandpassed audio
        if Ms >= .95
            disp("WARNING: Audio saturated")
            pause()
        end
        normalx = x_freq/Ms(1); %normalize bandpassed audio %needed?
        P_freq = bandpass(P,freq_bins,Fs); %bandpass signal (echosounder ping)
        MP_freq = max(P_freq); %maximum from template echosounder ping
        normalP_freq = P_freq/MP_freq; %normalize signal
        
        %%% attempt at moving search window loop %%%
         
        steps = [1:ceil((length(normalx)/N2))]'; %define number of steps using number of samples in audio and user defined width of search window
        
        b_freq = conj(normalP_freq(end:-1:1)); %calculate inverse conjugate of bandpassed signal
        fl1 = 60;
        y1 = filter(b_freq,1,normalx); %match filter
        
        
        [upall,loall] = envelope(y1,fl1,'analytic');
        [Peak_val_all, peak_loc_all] = findpeaks(upall,'MinPeakDistance',N);
        
        
        for n = 1:length(peak_loc_all)
            d_ping = peak_loc_all(n);
            %nstart=(n-1)*N2+1;
            %nend=nstart+N-1;
            %if nend > length(x_freq)
            %    nend = length(x_freq);
            %end      
            %window = nstart:nend;
            %win_x = x(window); %select out window
            %max_win = max(win_x);
            %norm_win_x = win_x/max_win; %normalize window
            %figure(5)
            %subplot(2,1,1)
            %plot(win_x)
            %y2 = filter(b_freq,1,win_x); %match filter
            %subplot(2,1,2)
            %plot(y2)
            %hold on
            %fl1 = 60;
            %[up1,lo1] = envelope(y2,fl1,'analytic');
            %M_up1 = max(up1);
            %norm_up1 = up1/M_up1;
            %[Peak_val, peak_loc] = findpeaks(up1,'MinPeakDistance',length(window)-2);
            %plot(up1)
            %hold on 
            %if ~isempty(peak_loc)
            %   plot(peak_loc,max(Peak_val),'r^')
               noise_window = d_ping - 1023:d_ping - 512;
               ping_window = d_ping - 255:d_ping + 256;
               if ping_window(end) > length(x_freq)
                   continue
               end 
               if ping_window(1) < 1
                  continue
               end
               if noise_window(end) > length(x_freq)
                   continue
               end 
               if noise_window(1) < 1
                  continue
               end
               ping = x_freq(ping_window);
               noise_sample = x_freq(noise_window);
               SNR = rssq(ping(:))/rssq(noise_sample(:));
               if SNR >= 0.5
               %[~,Max_ping] = max(Peak_val); %takes max ping
               ping_index = d_ping; %finds index in data
               ping_indices = [ping_indices; ping_index];
               ping_indices = unique(ping_indices);
               hold off
               end
               %pause() 
            %else
            %   hold off
            %   continue
         end
        %end                     % end search window loop
        
                                        %clean up plots for next loop
      %figure(1); clf;
      %figure(2); clf;
      slashIdx = strfind(PATH2WAV, '\'); 
      pathSegment = PATH2WAV(slashIdx(end)+1:end);
      file_n = split(pathSegment,'.');
      new_file = ['match_' num2str(frequencies(freq)) '_' char(file_n(1)) '.' char(file_n(2)) '.mat'];
      Path2Output = [PATH2WAV(1:slashIdx(end-1)) 'OUTPUT\'];
      if ~exist(Path2Output, 'dir')  
      mkdir(Path2Output);
      end
      save_path = string([Path2Output new_file]);
      %save(save_path, 'ping_indices')     
   end                         %end freq loop
end                             %end filelist loop

figure(8)
plot(x_freq)
hold on
plot(ping_indices,max(x_freq)*1.05,'r^')

% click length ~0.002s


