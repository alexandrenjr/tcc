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

SIGMA = -4/100;
TS = 0.5;
WN = 2;
SYSC = ss(Ac,Bc,Cc,Dc);

NY = 2*pi/(WN*TS);

if NY <= 4.86
  error('Curva não convexa!');
end

SYS = c2d(SYSC,TS,'tustin');

A = SYS.A;
B = SYS.B;
C = SYS.C;
D = SYS.D;

[n,m] = size(B);

P = sdpvar(n,n);
Z = sdpvar(m,n,'f');
options = sdpsettings('verbose',0);

F = [];
F = [F, (P>=0):'Positividade']; %#ok<*BDSCA>
F = [F, (lmiestabilidade(SIGMA,SYS,TS,P,Z)<=0):'Estabilidade relativa'];

k = 0;
No = pontoplanoz(0,WN,TS);
Ni = pontoplanoz(1,WN,TS);

pontos3 = [0 1];
pontos4 = pontos3;
vec3 = [Ni No];
vec4 = vec3;
theta = acos(abs(real(vec4(m+1))-Ni)/abs(vec4(m+1)-Ni));

Satual = 0.5*abs(Ni-real(No))*abs(No-real(No));

F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Ni)*P>=0):'Limitação à direita'];
F = [F, (lmisetorconico(Ni,theta,SYS,P,Z,'D')<=0):'Setor cônico direito'];

optimize(F,trace(P),options);

[primalres,dualres] = check(F);
primalres = sort(primalres,'ascend');
dualres = sort(dualres,'ascend');

while primalres(1) < 0 || dualres(1) < 0
  if k < length(vec3)-1
    k = k+1;
  else
    k = 1;
    vec3 = vec4;
    pontos3 = pontos4;
    clf
    hold on
    axis equal
  end
  hold on
  F = [];
  F = [F, (P>=0):'Positividade'];
  F = [F, (lmiestabilidade(SIGMA,SYS,TS,P,Z)<=0):'Establidade relativa'];

  zeta0 = pontos3(k);
  zeta1 = pontos3(k+1);
  pontos4 = sort([pontos4 (zeta0+zeta1)/2],'ascend');
  Vnew2 = pontoplanoz((zeta0+zeta1)/2,WN,TS);
  vec4 = sort([vec4 Vnew2],'ascend');

  F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Ni)*P>=0):'Limitação à direita'];

  for m=1:length(vec3)-1
    u2 = loc(vec4(m),vec4(m+1));
    theta = acos(abs(real(vec4(m+1))-u2)/abs(vec4(m+1)-u2));
    F = [F, (lmisetorconico(u2,theta,SYS,P,Z,'D')<=0):['Setor cônico direito ' num2str(m)]];
%     plot([real(vec4(m+1)),real(vec4(m)),real(u2),real(vec4(m+1))], ...
%       [imag(vec4(m+1)),imag(vec4(m)),imag(u2),-imag(vec4(m+1))],'m')
  end

  optimize(F,trace(P),options);
  [primalres,dualres] = check(F);
  primalres = sort(primalres,'ascend');
  dualres = sort(dualres,'ascend');
  
  Sant = Satual;
  Satual = polyshape(real(vec4),imag(vec4)).area;
  
  K = value(Z)/value(P);
  syscomp = ss(A+B*K,B,C+D*K,D,TS);
  pzmap(syscomp)

  if (Sant/Satual) < 1 && (Sant/Satual > 0.999999)
    efactivel = 0;
    disp(['Infactível por não atender o requisito wn=' num2str(WN)]);
    break;
  end
end

hold on
plot(real(taxadedecaimento(SIGMA,TS)),imag(taxadedecaimento(SIGMA,TS)),'m')
plot(real(vec4),imag(vec4),real(vec4),-imag(vec4))
K = value(Z)/value(P);
syscomp = ss(A+B*K,B,C+D*K,D,TS);
pzmap(syscomp,'r')
zgrid(-1,WN,TS)