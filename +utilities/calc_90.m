function xFilt = noDelayFilt(dFilter, x)
% Filter a signal using a digitalFilter object from the Signal Processing
% Toolbox, compensating for group delay introduced by the filter. This 
% function only works if the delay is not frequency-dependent (usually the
% case with FIR filters). It will generally NOT work with IIR filters like 
% the Butterworth filter.
%
% Last updated by Wilfried Beslin
% 2024-02-27