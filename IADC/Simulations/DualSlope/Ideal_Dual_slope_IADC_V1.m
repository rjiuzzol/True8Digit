%
%
% True8Digit dual-slope IADC simulation script
%

close all
clear all


%% internal script control variables
iFig=1; %figure counter (internal control)
dbg=0; %set Debugging variable
plt=1; %set plotting

%% Input parameters
% Integrator (based on HP3458 components)
R1=50e3;   % R= 50k
C1=330e-12; % C= 100 pF
R2=1e9; %switch off finite resistance
Vref=10; % 10 V voltage reference

% Sampling Parameters

Fclk=20e6; % system clock 20 MHz
Tclk=1/Fclk; % system clock period

Fint=10*1000; % Timer frequency for integrator circuit
Tint=1/Fint; % switch between Vin and Vref

Taper=82.5e-6;%5*R1*C1;%2e-6; % Aperture time
%Taper=5e-6;% As in LTSpice
%Taper=Tclk; %Taper minimun

if (Taper >= Tint)
  disp("Taper to high, please decrease it");
  msg=sprintf("TIMER: %f", Tint);
  disp(msg)
  msg=sprintf("Taper: %f", Taper);
  disp(msg)
  disp("Set Taper to a small value than Tint (TIMER)")
  exit;
endif


% input signal to be sampled
freq=1000; % input signal frequency
A0=1.00000;       % input signal amplitude
ph=0;      % input signal phase
Voff=0;   %DC offset voltage
pol=-1;   %signal polarity

fun = @(x)(A0*sin(2*pi*freq*x+ph)+Voff);
fref = @(x)Vref;


% Simulation execution time in input signal periods
nPeriods=1;
if (freq==0)
  nSamples=20;
  t_vec_end=1e-3;
else
  nSamples=Fint/freq;
  t_vec_end=nPeriods/freq;
endif

M=nPeriods*nSamples; %index


t = 0:1/Fint:t_vec_end; t=vec(t(1:end-1));
tclk=0:1/Fclk:t_vec_end; tclk=vec(tclk(1:end-1));
TIMER_t=0:Tint:t_vec_end; TIMER_t=vec(TIMER_t(1:end-1));
%swt = pulstran(tclk,TIMER_t+Taper/2,"rectpuls",Taper);

%% Sampling clock signal
arbclk=0:1/(2*Fclk):t_vec_end; arbclk=arbclk(1:end-1);
clksg=pulstran(arbclk,tclk,"rectpuls",Tclk/2);

% Signals generation
x1 = (pol)*fun(tclk);
x1=vec(x1);
xVref = ones(length(x1),1).*(-pol)*Vref;

% Auxiliaries Variables
TTs = [];


Tratio=round(Tint/Tclk);
TaperR=round(Taper/Tclk);



for k=0:M-1,
  idxi=k*Tratio+1;
  idxf=idxi+TaperR-1;

  TTs(k+1,3)=tclk(idxi);
  TTs(k+1,4)=tclk(idxf);
  TTs(k+1,5)=idxi;
  TTs(k+1,6)=idxf;

end

%% Solving using control toolbox functions
%int_response_ex;
% get system synthesis - WE NEED CONTROL TOOLBOX!
num=-1;
den=[R1*C1 R1/R2];
sys=getINTSS(num,den,dbg);

%% Solving in time-domain
%%
% Initial conditions
x0=0;%Voff/sys.c;

%%%%
%% Integration process
%%
%%
% Main LOOP
N=length(tclk);
sgpol=ones(length(TTs),1); %To keep signal polarity for later analysis

for ii=1:length(TTs),

  % Phase 1: input signal integration

  [xout(ii,:), xtout(ii,:), _nclkV(ii,:),tnc(ii,:)]=simulateINT(sys, x1(TTs(ii,5):TTs(ii,6)), ...
                                          tclk(TTs(ii,5):TTs(ii,6)), clksg(TTs(ii,5):TTs(ii,6)),...
                                          x0, [-10;10],'up');


  %check polarity
  sgpol(ii) = sgpol(ii)*sign(xout(end)./sys.c);

  %Phase 2: Voltage reference integration
  if ii < length(TTs)
    [vout(ii,:), vtout(ii,:), _nclkR(ii,:),tnc(ii,:)]=simulateINT(sys, xVref(TTs(ii,6)+1:TTs(ii+1,5)-1),...
                                            tclk(TTs(ii,6)+1:TTs(ii+1,5)-1), clksg(TTs(ii,6)+1:TTs(ii+1,5)-1), ...
                                            abs(xout(end)./sys.c), [0;10],'down');

  else
    [vout(ii,:), vtout(ii,:), _nclkR(ii,:),tnc(ii,:)]=simulateINT(sys, xVref(TTs(ii,6)+1:end), tclk(TTs(ii,6)+1:end), ...
                                            clksg(TTs(ii,6)+1:length(x1)),abs(xout(end)./sys.c), [0;10], 'down');

  endif


end


%% Reshaping
xx=reshape(xout',1,((M)*TaperR));
xt=reshape(xtout',1,((M)*TaperR));
nClkV=reshape(_nclkV', 1, length(_nclkV));
nClkV = vec(nClkV);

xv=reshape(vout',1,(M*(Tratio-TaperR)));
xtv=reshape(vtout',1,(M*(Tratio-TaperR)));
nClkR=reshape(_nclkR', 1, length(_nclkR));
nClkR = vec(nClkR);


% Counter ratio, signal reconstruction
countR=[nClkR nClkV nClkR./nClkV (nClkR./nClkV).*Vref.*(sgpol)];

%%
% output signal

sgout=countR(:,4);

% Taper correction
%Kaper=sin(pi*freq*Taper)/(pi*freq*Taper);
%sgout_cor=sgout/Kaper;

%% Analysis in the frequency domain
% Fourier transform
%
f_vec=[0:1/M:1].*Fint;
f_vec=vec(f_vec);
X2=fft(sgout);


Mag=2*abs(X2)/length(X2);
Phase=angle(X2).*(180/pi);

fbin=nPeriods+1;
As=Mag(fbin);
Phi=Phase(fbin);


%%
% Display results
msg=sprintf('Sampling Rate (Hz): %f', Fint);
disp(msg);
msg=sprintf('Aperture Time (s): %e', Taper);
disp(msg);
msg=sprintf('Input signal parameters');
disp(msg);
msg=sprintf('Amplitude (V): %f',A0);
disp(msg);
msg=sprintf('Offset (V): %f',Voff);
disp(msg);
msg=sprintf('Frequency (Hz): %f',freq);
disp(msg);

disp('Results')
disp('Time domain')
msg=sprintf('Sampled Amplitude (V): %f',std(sgout,1)*sqrt(2));
disp(msg);
disp('Frequency domain')
msg=sprintf('Sampled Amplitude (V): %f',As);
disp(msg);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ploting


if plt

  figure(iFig++,"position",[0*640 300 640 480]);
  plot(xtv,xv,'bo');
  hold on;
  plot(xt,xx,'r*');grid
  hold off;
  title('Integration Process')
  legend('Voltage Reference','Input signal');

  figure(iFig++,"position",[1*640 300 640 480]);
  plot(t(1:end),sgout,'r*');grid;
  title('Sampled signal')
  xlabel('time (s)')
  ylabel('Amplitude (V)')

  figure(iFig++,"position",[2*640 300 640 480]);
  plot(f_vec(1:end/2),20*log10(Mag(1:end/2)), 'r*-');
  grid;
  xlabel('Frequency (Hz)')
  ylabel('Amplitude (V)')
  title('Sampled signal spectrum')


end



