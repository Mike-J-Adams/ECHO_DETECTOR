function y=Pascals2wav(p)
%
% see: ACADIA/RICHARD/Matlab_For_Brian/RANGE/SOUND/PORPOISE/OCEAN_SONICS/REPORTS/manual.dvi
%
% y=Count/2^(bits-1)
%
% example:
%  [yy,fs,nbits]=wavread('M05yy161555.wav');
%  p=wav2Pascals(yy);
% load tSS.mat
% SSP=wav2Pascals(SS);
% SS=Pascals2wav(SSP);

VpK=3;
y=p/(VpK*10^(169/20)*10^-6);   % Pascals
