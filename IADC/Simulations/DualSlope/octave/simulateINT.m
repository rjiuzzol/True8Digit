function [y, t, nClk, tClki]=simulateINT(sys, u, tclk, clksg, x0, x0p, intphase)

[urows, ucols] = size (u);

if nargin==7 && strcmp(intphase,'up')
  ph=0;% running up
elseif nargin==7 && strcmp(intphase,'down')
  ph=1; %running down
endif


%sys1 = ss(A,B,C,D);
method = "foh";
dt = abs (tclk(end) - tclk(1)) / (urows - 1);   # assume that t is regularly spaced
%t = vec (linspace (tclk(1), tclk(end), urows));
t = vec (tclk);
%dt = abs (t1(end) - t1(1)) / (urows - 1);   # assume that t is regularly spaced

sysc = c2d (ss (sys), dt, method);           # convert to discrete-time model (in ss for accuracy)
[A, B, C, D] = ssdata (sysc);
[p, m] = size (D);                            # number of outputs and inputs
n = rows (A);                                 # number of states
x = ones(1,1).*x0;
x = x - sysc.userdata * u(1,:)';

nClku=0; %number of Runup clock pulses
nClkd=0; %number of Rundown clock pulses
for i=1:urows
    y(i,:) = C*x + D*u(i,:).';
    x = A* x + B * u(i,:).';

    k = y(i,:) > x0p(2);
    if any(k)
      y(i,k) = x0p(2);
      if ph==0
        disp('Overloading');
        %break;
      endif
    endif
    k = y(i,:) < x0p(1);
    if any(k)
        y(i,k)= x0p(1);
        nClkd++;
      if ph==0
        disp('Overloading');
        %break;
      endif

    endif

    %if (y(i,:)<x0p)
    %  y(i,:)=x0p;
    %  nClk--;
    %endif

    nClku++; % Number of clock cycles to integrate the Voltage
end
tClki=sum(clksg==1);
%if nClk < 0
tClki = tClki - (nClkd/2);
%endif


nClk = nClku-nClkd;

endfunction

