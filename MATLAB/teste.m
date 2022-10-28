clear;
clc;

zeta = 0.9;
Ts = 1;
a = double(realwn(zeta,Ts))/(Ts*sqrt(1-zeta^2));
b = pontoplanoz(zeta,a,Ts);
hold
plot([real(b),real(b)],[0,imag(b)])
zgrid(zeta,-1,-1)

% theta = 0:1e-6:pi;
% hold
% plot(theta,cos(theta))
% plot(theta,-exp(zeta*(theta-pi)/sqrt(1-zeta^2)))