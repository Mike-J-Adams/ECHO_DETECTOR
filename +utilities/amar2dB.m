function [p] = amar2dB(x,Gain,Hs)
% convert digital values from AMAR to linear pressure units 
% 
% Input:
% x = digital AMAR values
% Gain = ADC gain for G4
% Hs = calculated hydrophone sensitivity from pistonphone calibration
% 
% Output:
% p = data in calibrated linear pressure units
% pdB = data in calibrated SPL in dB
%
% Last updated by Mike Adams
% 2024-05-28

H = 10^(Hs/20);% convert sensitivity to linear
M = 5.36*10^-7; % conversion factor

p = ((double(x)*M)./(Gain*H));
