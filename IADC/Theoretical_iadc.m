
%%
% Sampling a Sinus with a Theoretical IADC based on HP3458
%
% Authors: Ricardo Iuzzolino
%
% Date: 1-February-2024
%
% version: 1.0


close all;
clear all;

% Input signal parameters
freq=62.5;
nPeriod=8;
nMuestras=16;
A0=0.5*8;

% Sampling parameters (calculated according the input signal to have a coherent sampling)
fs = nMuestras*freq;
t = 0:1/fs:nPeriod/freq; t=t(1:end-1);
ts=0:1/(nMuestras*fs):nPeriod/freq; ts=ts(1:end-1);

% Input signal
x1 = A0*sin(2*pi*freq*t);

% Sampling time
Ts = 1/fs;

% Aperute time (HP3458)
SAFETIME = 25e-6;
Ta = Ts - SAFETIME;


%% Input signal to be sampled (edit the fun(x))
fun = @(x)A0*sin(2*pi*freq*x);
q=zeros(nMuestras*nPeriod,1);

for k=0:nPeriod*nMuestras-1,
q(k+1) = integral(fun,k*Ts,k*Ts+Ta)/Ta;
end

%
% Fourier transform
%
X1=fft(x1);
X2=fft(q);


Mag=2*abs(X2)/length(X2);
Phase=angle(X2).*(180/pi);

disp('Results')
fbin=nPeriod+1;
As=Mag(fbin)
Phi=Phase(fbin)


% plots
figure
plot(t, q,'r*-');
xlabel('time (s)')
ylabel('amplitude (au)')
title('sampled signal by IADC')

figure
subplot(211),plot((Mag(1:end/2)), 'r*-');
subplot(212),plot(Phase(1:end/2),'ro-');

KSinc=sin(pi*freq*Ta)/(pi*freq*Ta)
KPhase = angle(exp(-j*pi*freq*Ta))*180/pi


%% Corrected results for t_aper

As/KSinc/sqrt(2)
Phi+KPhase
