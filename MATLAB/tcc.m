clear;
clc;

R1 = 1;
C1 = 5;
R2 = 1;
C2 = 5;
Ceq = (C1*C2)/(C1+C2);

Ac = [-inv(R1*Ceq) -inv(R1*C2); inv(R2*C2) -inv(R2*C2)];
Bc = [inv(R1*C1); 0];
Cc = [R1 0];
Dc = zeros;

SIGMA = -4/80;
ZETA = 0.6;
TS = 1.9;
WS = (2*pi)/TS;
WN = 0.5;
NY = WS/WN;

SYSC = ss(Ac,Bc,Cc,Dc);
SYS = c2d(SYSC,TS,'tustin');

K = factibilidade(SYS,TS,SIGMA,ZETA,WN,'P',1);

SYSCOMP = ss(SYS.A+SYS.B*K,SYS.B,SYS.C+SYS.D*K,SYS.D,TS);
figure
step(SYSCOMP)
