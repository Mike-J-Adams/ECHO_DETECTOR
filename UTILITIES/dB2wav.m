function y=dB2wav(L)
%
% see: ACADIA/RICHARD/Matlab_For_Brian/RANGE/SOUND/PORPOISE/OCEAN_SONICS/REPORTS/manual.dvi
% y=Count/2^(bits-1)

VpK=3;

% L=20*log10(VpK*y)+169;   % L dB    % formula for wav2dB

y=(L-169)/20;

y=10.^y;

y=y/VpK;
