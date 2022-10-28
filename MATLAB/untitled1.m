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

TS = 1;
SIGMA = -4/5;
ZETA = 0.6;
wnv = 0:1e-2:pi/(TS*sqrt(1-ZETA^2));

sysc = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sysc,TS,'tustin');

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

[n,m] = size(B);

P = sdpvar(n,n);
Z = sdpvar(m,n,'f');
options = sdpsettings('verbose',0);

F = [];
F = [F, (P>=0):'Positividade']; %#ok<*BDSCA>
F = [F, (lmiestabilidade(SIGMA,sysd,TS,P,Z)<=0):'Estabilidade relativa'];

k = 0;
Vo = pontoplanoz(ZETA,0,TS);
Vf = pontoplanoz(ZETA,pi/(sqrt(1-ZETA^2)*TS),TS);
Vfbarra = pontoplanoz(ZETA,double(realwn(ZETA,TS)),TS);

pontos1 = [0 pi/(1.5*sqrt(1-ZETA^2)*TS)];
pontos2 = pontos1;
vec1 = [Vo Vfbarra];
vec2 = vec1;

Satual = polyshape(real([vec2 Vf]),imag([vec2 Vf])).area;

F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Vf)*P>=0):'Limitação à direita'];
F = [F, (lmisetorconico(loc(Vo,Vfbarra),acos(abs(Vo-real(Vfbarra))/abs(Vo-Vfbarra)), ...
  sysd,P,Z,'E')<=0):'Setor cônico Esquerdo'];

sol = optimize(F,trace(P),options);

[primalres,dualres] = check(F);
primalres = sort(primalres,'ascend');
dualres = sort(dualres,'ascend');

plot(real(vec2),imag(vec2))

while primalres(1) < 0 || dualres(1) < 0
  if k < length(pontos1)-1
    k = k+1;
  else
    k = 1;
    vec1 = vec2;
    pontos1 = pontos2;
  end
  
  F = [];
  F = [F, (P>=0):'Positividade'];
  F = [F, (lmiestabilidade(SIGMA,sysd,TS,P,Z)<=0):'Estabilidade relativa'];

  Vnew1 = pontoplanoz(ZETA,(pontos1(k)+pontos1(k+1))/2,TS);
  pontos2 = sort([pontos2 (pontos1(k)+pontos1(k+1))/2],'descend');
  vec2 = sort([vec2 Vnew1],'descend');

  F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Vf)*P>=0):'Limitação à direita'];
  
  for m=1:length(vec1)-1
    u1 = loc(vec2(m),vec2(m+1));
    if u1 < 0
      phi = acos(abs(real(vec2(m+1))-u1)/abs(vec2(m+1)-u1));
      F = [F, (lmisetorconico(u1,phi,sysd,P,Z,'D')<=0):['Setor cônico direito ' num2str(k)]];
    else
      phi = acos(abs(u1-real(vec2(m+1)))/abs(u1-vec2(m+1)));
      F = [F, (lmisetorconico(u1,phi,sysd,P,Z,'E')<=0):['Setor cônico esquerdo ' num2str(k)]];
    end
  end

  sol = optimize(F,trace(P),options);
  [primalres,dualres] = check(F);
  primalres = sort(primalres,'ascend');
  dualres = sort(dualres,'ascend');
  
  Sant = Satual;
  Satual = polyshape(real(vec2),imag(vec2)).area;

  if (Sant/Satual) < 1 && (Sant/Satual > 0.999999)
    disp('Infactível!');
    break;
  end
  hold on
  pgon = polyshape([real([real(Vf) vec2])],[imag([0 vec2])]);
  plot(pgon);
end

hold on
axis equal
xline(real(Vf),'m')
plot(real(taxadedecaimento(SIGMA,TS)),imag(taxadedecaimento(SIGMA,TS)),'m')
plot(real(vec2),imag(vec2),real(vec2),-imag(vec2))
K = value(Z)/value(P);
syscomp = ss(A+B*K,B,C+D*K,D,TS);
pzmap(syscomp,'r')
zgrid(ZETA,-1,-1)