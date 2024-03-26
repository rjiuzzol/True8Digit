Ideal Dual-Slope ADC simulation files
% Sampling Parameters
NCounts=120000000; % 8 1/2 digits
Fclk=20e6; % system clock 20 MHz
Tclk=1/Fclk; % system clock period

Fint=10*1000; % Timer frequency for integrator circuit
Tint=1/Fint; % switch between Vin and Vref

Taper=5*R1*C1;%2e-6; % Aperture time
