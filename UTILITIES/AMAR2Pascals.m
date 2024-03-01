function p = AMAR2Pascals(x,V,S)
%  𝑝[μPa]=DV[bits]×𝑀[V/bit]𝐺[V/V]×𝐻[V/μPa]
%   convert AMAR voltage to pascals (uPa)

%x=10000;
%S=-164;

M = 5.36*10^-7; %conversion factor for AMAR G4

if exist('V')==0
    V = 5;  %factory default voltage G4
end
if exist('S')==0
    S = -200; %Geospectrum sensitivity M36-100
end

S = 10^(S/20);

p = (x*M)/(V*S);
end

