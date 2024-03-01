function [ishift,Cenv,max_env,y]=find_peak_Hilbert(x,X,A,n1,n2,npeak,swplt)
%
% Calculates the envelope using a Hilbert transform. Finds where its
% peak is and correlates the envelope with ClickAmplitude.
%
% x is the signal that has been band-pass filtered to pass porpoise clicks
% X is the FFT of x
% A is the ClickAmplitude function which is maximum at imaxA=65
% length(x)=length(A)=NFFT (=128)
% The algorithm considers A(n1:n2) where n2>n1>=31 and n2<=NFFT-30
% swplt    1 for plots, 0 to turn plots off
%
% ishift   number of indices that x must be shifted in order to match A
% Cenv     correlation between A and the amplitude of shifted x
% max_env  maximum value of the click envelope
% y=x(n1+ishift:n2+ishift);   % select the shifted segment of data
%
global Octave

NFFT=length(x);

[F,E]=log2(NFFT); if F~=0.5, error('NFFT must be a power of 2'), end

if length(X)~=NFFT | length(A)~=NFFT
 error('dimensions of X or A are wrong in find_peak_Hilbert')
end

% reconstruct the signal with a pi/2 phase shift
Y=X; Y(2:NFFT/2)=-i*Y(2:NFFT/2); Y(NFFT/2+2:end)=i*Y(NFFT/2+2:end);
Y(1)=0; Y(NFFT/2+1)=0;
y=ifft(Y,NFFT);
y=real(y);      % Octave patch

% calculate the envelope
env = sqrt(x.^2+y.^2);

%disp(num2str([[0:127]' x y env]))
%figure(1),plot([1:128],x,'r',[1:128],y,'b')

%[max_env,ipeak]=max(env);   
%ishift=ipeak-65;
%disp(['ipeak=',num2str(ipeak)])

% smooth the envelope
%env(n1:n2)=(env(n1-1:n2-1)+2*env(n1:n2)+env(n1+1:n2+1))/4;
env(2:end-1)=(env(1:end-2)+2*env(2:end-1)+env(3:end))/4;

[max_env,ipeak]=max(env(n1:n2));   

ishift=ipeak+n1-1 -npeak;

ishift=min(max(ishift,-n1+1),NFFT-n2);

if swplt==1
 h=figure(3); set(h,'position',[60 527 560 420]),plot(x),xlabel('x(1:NFFT)','fontsize',14)
 title(['ipeak=',num2str(ipeak),' initial ishift = ',num2str(ishift)],'fontsize',14)
end
% disp(num2str([A(n1:n2) env(n1+ishift:n2+ishift)]))

%q=ver; 
%if strcmp(q(1).Name,'Matlab')
if Octave
 C=corr(A(n1:n2),env(n1+ishift:n2+ishift));       % Octave
else
 C=corrcoef(A(n1:n2),env(n1+ishift:n2+ishift));   % Matlab
end

% Matlab returns a matrix from corrcoef whereas Octave returns a scalar
%q=ver; if q(1).Name(1)=='M', Cenv=C(1,2); else, Cenv=C; end
if length(C(:))>1, Cenv=C(1,2); else, Cenv=C; end

y=x(n1+ishift:n2+ishift);   % select the shifted segment of data

if swplt==1
 figure(4),plot([n1:n2],max_env*A(n1:n2),'r',[n1:n2],env(n1+ishift:n2+ishift),'b-',[n1:n2],y,'g')
 legend('std Click','env','y','location','northwest')
 xlabel(['n1=',num2str(n1),'  n2=',num2str(n2),'  ishift=',num2str(ishift)],'fontsize',14)
 title(['Created by find\_peak\_hilbert.m Cenv=',num2str(Cenv)],'fontsize',14)
end
