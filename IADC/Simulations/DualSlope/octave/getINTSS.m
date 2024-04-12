function [sys] = getINTSS(num, den, dbg)

if dbg
  disp('*** getINTSS ')
endif

if ~exist('tf2ss')
  pkg load control;
endif


H=tf(num,den);
[A,B,C,D]=tf2ss(num,den);
x0 = zeros(1,1);

[int.zeros,int.poles,int.k] = zpkdata(H,"v");

sys = ss(A,B,C,D);

%%%% Debug control
if dbg

  disp('Transfer function')
  H

  disp('ZPK variables')
  int

  disp('State variables')
  A
  B
  C
  D

  disp('SYS control system')
  sys
endif

endfunction


