function L=wav2dB(y)
%
% see: ACADIA/RICHARD/Matlab_For_Brian/RANGE/SOUND/PORPOISE/OCEAN_SONICS/REPORTS/manual.dvi
% y=Count/2^(bits-1)

VpK=3;
L=20*log10(VpK*y)+169;   % L dB

