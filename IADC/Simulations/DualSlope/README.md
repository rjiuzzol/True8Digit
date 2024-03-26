Ideal Dual-Slope ADC simulation files

% Octave
To simulate using GNU Octave

first load signal and control packages using the command line:

pkg load signal

pkg load control


Then open the file "Ideal_Dual_Slope.m" with the built-in editor and configure the signal and sampling parameters as follow:

% Sampling Parameters

Fclk=20e6; % system clock 20 MHz

Fint=10*1000; % Timer frequency for integrator circuit

Taper=5*R1*C1;%2e-6; % Aperture time


Finally set the input signal parameters:

% input signal to be sampled

freq=1000; % input signal frequency

A0=1.0000054;       % input signal amplitude

ph=0;      % input signal phase

Voff=A0;   %DC offset voltage

pol=-1;   %signal polarity


