clear;
clc;

R1 = 1;
C1 = 5;
R2 = 1;
C2 = 5;
Ceq = (C1*C2)/(C1+C2);

A = [-inv(R1*Ceq) -inv(R1*C2); inv(R2*C2) -inv(R2*C2)];
B = [inv(R1*C1); 0];
C = [R1 0];
D = zeros;
Bd = [.1; 0];
Dd = zeros;

sigma = -4/100;
zeta = 0.7;
Ts = 2.9;
wn = 2;
Ny = 2*pi/(wn*Ts);

sysc = ss(A,B,C,D);
sysd = c2d(sysc,Ts,'tustin');

K = factibilidade(sysd,Ts,sigma,zeta,wn,'C',true);

% syscomp = ss(sysd.A+sysd.B*K,sysd.B,sysd.C+sysd.D*K,sysd.D,Ts);
% 
% hold on
% axis equal
% pzmap(syscomp,'r')
% 
% figure
% step(syscomp,'r')