clear;
clc;

C1 = 5;
C2 = 5;
R1 = 1;
R2 = 1;

Ceq = (C1*C2)/(C1+C2);

Ac = [-inv(R1*Ceq) -inv(R1*C2); inv(R2*C2) -inv(R2*C2)];
Bc = [inv(R1*C1); 0];
Cc = [R1 0];
Dc = zeros;

ts = 50;
SIGMA = -4/ts;
ZETA = 0.5;
TS = 10;
WS = (2*pi)/TS;
WN = 0.1;
NY = WS/WN;

SYSC = ss(Ac,Bc,Cc,Dc);
SYS = c2d(SYSC,TS,'tustin');

K = factibilidade(SYS,TS,SIGMA,ZETA,WN,'P');

figure
SYSCOMP = ss(SYS.A+SYS.B*K,SYS.B,SYS.C+SYS.D*K,SYS.D,TS);
impulse(SYS,SYSCOMP)
