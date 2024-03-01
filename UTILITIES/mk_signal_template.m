%mk_signal_template.m

% read in strong pings and generate template

clear
close all

Fs = 256000; %sample rate

[file,path] = uigetfile('*.wav'); %User selects file clip containing strong signal
PATH2PING = [path file];
SP = audioread(PATH2PING);
MSP = max(SP);
normalSP = SP/MSP(1); %normalize to max
fl1 = 60;
[up,lo] = envelope(normalSP,fl1,'analytic');
figure(1)
subplot(4,1,1)
plot(normalSP)
hold on
plot(up)

% find local peaks
[Peak_val_main, peak_loc_main] = findpeaks(up,'MinPeakProminence',0.95);
hold on 
if ~isempty(peak_loc_main)
    plot(peak_loc_main,max(Peak_val_main)*1.05,'r^')
end

hold off

%extract main signal
window_size = 256;
main_sig_index = peak_loc_main-window_size:peak_loc_main+window_size-1;
main_sig = SP(main_sig_index);
subplot(4,1,2)
plot(main_sig)
%fl1_main = 40;
%[Main_up,Main_lo] = envelope(main_sig,fl1_main,'analytic');
Main_env_peak = envelope(main_sig,10,'peak'); % peak envelope might be better for this extent
hold on
%plot(Main_up)
plot(Main_env_peak)
hold off

Main_peak_X = fft(main_sig);
%figure(3)
%plot(Main_peak_X)
subplot(4,1,3)
pspectrum(main_sig,Fs)
subplot(4,1,4)
pspectrum(main_sig,Fs,'spectrogram')

%remove main signal to get reflection
rm_main = SP;
rm_main(main_sig_index) = 0;
rm_main_max = max(rm_main);
normal_rm_main = rm_main/rm_main_max(1); %normalize to max
flR = 60;
[upR,loR] = envelope(normal_rm_main,flR,'analytic');
[Peak_val_R, peak_loc_R] = findpeaks(upR,'MinPeakProminence',0.95);
window_mod = 2; %widen window
R_peak_index = peak_loc_R(1)-window_size*window_mod:peak_loc_R(1)+window_size*window_mod-1;
Ref_sig = rm_main(R_peak_index);
Ref_env_peak = envelope(Ref_sig,10,'peak'); % peak envelope might be better for this extent
figure(2)
subplot(4,1,1)
plot(normal_rm_main)
hold on
if ~isempty(peak_loc_main)
    plot(peak_loc_R,max(Peak_val_R)*1.05,'r^')
end
subplot(4,1,2)
plot(Ref_sig)
hold on
plot(Ref_env_peak)
hold off
subplot(4,1,3)
pspectrum(Ref_sig,Fs)
subplot(4,1,4)
pspectrum(Ref_sig,Fs,'spectrogram')

%audiowrite("INPUT\AMAR538_STRONG_PING_TEMPLATE_1.wav",main_sig, Fs,'BitsPerSample',24)
%audiowrite("INPUT\AMAR538_REFLECTED_PING_TEMPLATE_1.wav",Ref_sig, Fs,'BitsPerSample',24)
