%Create_template.m

clear
close all

%Fs = 256000; ek60
Fs = 8000;

%N=1024;                 % sample length

%N=Fs;   %N = Fs pull of 1 seconds of data
%N2=N/2;                 % step size = 50% overlap .5 seconds
%frequencies =[18000, 38000, 50000, 70000, 120000]; %echosounder EK60
frequencies =[20 300]; %Minke start stop frequencies
%bandpass_width = 5000; % +- width of bandpass filter %Ek60
bandpass_width = []; % +- width of bandpass filter %Minke
SNR_threshold = 2;
fl1= 60;

%set up rough window
%window = 0.004; %seconds %EK60
window = 0.25; %seconds %Minke
window_samps = window*Fs;

%Load INPUT pings 
PATH2INPUT = "E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\MINKE\";
input_list = dir(strcat(PATH2INPUT, "*.wav"));

PINGS = [];
SNR1 = [];
SNR2 = [];
RMS90 = [];
DUR90 = [];
SAMPS90 = [];

%load strong pings
for p = 1:length(input_list)
PATH2PING = [input_list(p).folder, '\', input_list(p).name];
[P] = audioread(PATH2PING);
[MP, MPI] = max(P);
absP = abs(P);
[AMP, AMPI] = max(absP);
%signal
PT_window = MPI - ceil(window_samps/2):MPI + ceil(window_samps/2)-1;
PW = P(PT_window);
figure(1)
subplot(2,1,1)
plot(PT_window,PW)
hold on
ylim([-1 1])
cumenergy = [0];
absPW = abs(PW);
cumPW = cumsum(absPW);

%%%%% think about this whole section more... cumulative energy is better.
%%%%% RMS is meaningless
maxPW = max(cumPW); %117
minPW = min(cumPW);
PW951 = interp1(cumPW, cumPW, maxPW*0.95,'previous'); 
PW952 = interp1(cumPW, cumPW, maxPW*0.95,'nearest');  
if PW951 <= maxPW*0.95
    PW95 = PW952;
elseif PW951 > maxPW*0.95
    PW95 = PW951;
end  
PW051 = interp1(cumPW, cumPW, maxPW*0.05,'next'); 
PW052 = interp1(cumPW, cumPW, maxPW*0.05,'nearest');  
if PW051 >= maxPW*0.05
    PW05 = PW052;
elseif PW051 < maxPW*0.05
    PW05 = PW052;
end
%%%%%

index95 = find(cumPW == PW95);
index05 = find(cumPW == PW05);
xline(index05+min(PT_window))
xline(index95+min(PT_window))
hold off
PW90_window = [index05+min(PT_window):index95+min(PT_window)];
PW90 = P(PW90_window);
rms90 = rms(PW90);
RMS90 = [RMS90, rms90];
dur90 = length(PW90)/Fs;
DUR90 = [DUR90, dur90];
samps90 = length(PW90);
SAMPS90 = [SAMPS90, samps90];
%for s = 1:length(PW)
%    temp = PW(s);
%    temprms = rms(temp);
%    if s==1
%        cumRMS = [rms];
%    else
%        cumrms = cumRMS(s-1)+rms;
%        cumRMS = [cumRMS; cumrms];
%    end
%nClickMax = round(355*(Fs/250000));
%dur = length(PW90)/Fs*1000;
%nClickActual = round((dur/1000)*Fs);
%nClick = min([nClickActual,nClickMax]);   
%rmsClickLin = sqrt(sum(PW90(101:101+nClick).^2)/nClick);
%rmsSignal = 20*log10(rmsClickLin);  
% noise
noise_window = PT_window-(window_samps/2+1024);
noise_sample = P(noise_window);
noise90_window = PW90_window - (window_samps/2+1024);
noise90 = P(noise90_window);
subplot(2,1,2)
plot(noise_window,noise_sample)
hold on
xline(max(noise90_window))
xline(min(noise90_window))
ylim([-1 1])
hold off
SN_window = [min(noise_window)-512:max(PT_window)+512];
SN = P(SN_window);
figure(2)
plot(SN_window,SN,PT_window,PW,'g',noise_window,noise_sample,'r')
ylim([-1 1])
%nNoise = numel(noise_sample);     
%rmsNoiseLin = sqrt(sum(noise_sample.^2)/nNoise);
%rmsNoise = 20*log10(rmsNoiseLin);
% SNR
%snr1 = rmsSignal - rmsNoise;
%SNR1 = [SNR1, snr1];
snr1 = snr(PW90,noise90);
SNR1 = [SNR1, snr1];
snr2 = snr(PW,noise_sample);
SNR2 = [SNR2, snr2];

%ping_envelope = envelope(PW90,300,'analytic');
%figure(1)
%subplot(2,1,1)
%hold on
%plot(PW90_window,ping_envelope,'r')
%hold off
PINGS = [PINGS, PW];

end

Max_samps90 = max(SAMPS90);
mean_ping = mean(PINGS,2);
figure(3)
plot(mean_ping)
PATH2MEANPING = strcat(PATH2INPUT, 'Mean_Ping\MEAN_MINKE.wav');




audiowrite(PATH2MEANPING, mean_ping, Fs,'BitsPerSample',24)
