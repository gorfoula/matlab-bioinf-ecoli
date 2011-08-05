t = 0:900; A = 1000; a = 0.005; b = 0.005;
z1 = A*exp(-a*t);
z2 = sin(b*t);
[haxes,hline1,hline2] = plotyy(t,z1,t,z2,'semilogy','plot');

[haxes,hline1,hline2] = plotyx(t,z1,t,z2,'semilogy','plot');

% semilogy/x