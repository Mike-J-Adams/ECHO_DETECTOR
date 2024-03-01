%env2sig.m

%testing taking the peak evnvelope of the mean bandpassed ek60 ping and generating single frequency sine wave ping 
clear 
close all

[P, Fs] = audioread('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\COC\Mean_Ping\MEAN_PING.wav'); %ek60 ping
%set up Window
window = 0.004; %seconds
window_samps = window*Fs;
frequencies=18000;
bandpass_width = 2500; % +- width of bandpass filter 

[MP, MPI] = max(P);
PT_window = MPI - ceil(window_samps/2):MPI + ceil(window_samps/2)-1;
PW = P(PT_window);
freq_bins = [frequencies-bandpass_width frequencies+bandpass_width];
P_freq = bandpass(PW,freq_bins,Fs); %bandpass signal (echosounder ping)
MP_freq = max(P_freq); %maximum from template echosounder ping
normalP_freq = P_freq/MP_freq; 
figure(1)
plot(normalP_freq)
hold on
[up, lo] = envelope(normalP_freq,300);
plot(up)
plot(lo)
hold off
%example code to generate sinwave
 % Sampling frequency (samples per second) 
 dt = 1/Fs; % seconds per sample 
 StopTime = length(up)/Fs; % seconds 
 t = (0:dt:StopTime)'; % seconds
 t(end,:) = []; 
 data = sin(2*pi*frequencies*t); 
 test = data.*up;
 figure(2)
 plot(t,data)
 figure(3)
 plot(test)
 
 %audiowrite('E:\BW_ECHO_EXPERIMENT\MATLAB\ECHO_DETECT\INPUT\18kHzPing.wav', test, Fs,'BitsPerSample',24)
 
