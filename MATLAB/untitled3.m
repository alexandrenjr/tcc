clear;
clc;

Ts = 1;
zetav = 0:0.01:1;
y = 0:0.1:0.5;
% epsilon = pi/(2*Ts):0.01:3/2*pi/Ts;

% hold on

%% S
% for sigma = -0.05:-0.1:-0.5
%   xline(sigma,'--k');
% end
% 
% wn = 0:0.01:0.5;
% 
% for zeta = 0.2:0.1:0.8
%   a = pontoplanos(zeta,wn);
%   plot(real(a),imag(a),'--k',real(a),-imag(a),'--k')
% end
% 
% zeta = 0:0.01:1;
% 
% for wn = 0.1:0.05:0.5
%   a = pontoplanos(zeta,wn);
%   plot(real(a),imag(a),'--k',real(a),-imag(a),'--k')
% end

%% Z

% for sigma = 0.2:0.1:1.2
%   a = taxadedecaimento(sigma,Ts);
%   plot(real(a),imag(a),'--k')
% end

zeta = 0.5;
wnv = 0:1e-3:pi/(Ts*sqrt(1-zeta^2));
a = pontoplanoz(zeta,wnv,Ts);
pgon = polyshape(real(a),imag(a));
disp(pgon.area);

% for zeta = 0.1:0.1:1
%   wnv = 0:pi/(Ts*sqrt(1-zeta^2)*100):pi/(Ts*sqrt(1-zeta^2));
%   a = pontoplanoz(zeta,wnv,Ts);
%   plot(real(a),imag(a),'--k',real(a),-imag(a),'--k')
% end

% for wn = 0:0.1*pi/Ts:pi/Ts
%   a = pontoplanoz(zetav,wn,Ts);
%   plot(real(a),imag(a),'--k',real(a),-imag(a),'--k')
% end

wn = 1;
Ny = 2*pi/(wn*Ts);
% 
% a = 1-exp((-2*pi)/Ny);
% b = (a^2*sin((-2*pi)/Ny))/(sqrt(a^2-(cos((-2*pi)/Ny)-1)^2));

% res = pontoplanoz(zetav,wn,Ts);
% plot(real(res),imag(res),'--k')
% 
% epsilon = 0:0.01:1;
% t = exp(-abs(0)*Ts)*exp(1i.*epsilon.*Ts);
% plot(real(t),imag(t),'--k')

% syms u v
% fimp = fimplicit((u-1)^2/a^2+(a^2-(cos((-2*pi)/Ny)-1)^2)*v^2/(a^2*sin((2*pi)/Ny)^2)==1, ...
%   [exp(-2*pi/Ny) real(pontoplanoz(0,wn,Ts)) -imag(pontoplanoz(0,wn,Ts)) imag(pontoplanoz(0,wn,Ts))], ...
%   'Color','k');

% T = [];
% T(:,1) = fimp.XData';
% T(:,2) = fimp.YData';
% writematrix(T, 'data.txt','Delimiter',' ');

% zgrid(-1,wn,1)
% matlab2tikz('G:\Meu Drive\TCC\tcc\Monografia\figuras\aproximacao_eliptica_wn.tikz')