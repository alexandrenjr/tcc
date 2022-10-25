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

Ts = 2.748180854900330;
sigma = -0.040000000000000;
zeta = 0.5;
wnv = 0:1e-2:pi/(Ts*sqrt(1-zeta^2));

sysc = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sysc,Ts,'tustin');

A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

[n,m] = size(B);

P = sdpvar(n,n);
Z = sdpvar(m,n,'f');
options = sdpsettings('verbose',0);

Fp = [];
Fp = [Fp, (P>=0):'Positividade']; %#ok<*BDSCA>
Fp = [Fp, (lmiestabilidade(sigma,sysd,Ts,P,Z)<=0):'Estabilidade relativa'];

k = 0;
Vo = pontoplanoz(zeta,0,Ts);
Vf = pontoplanoz(zeta,pi/(sqrt(1-zeta^2)*Ts),Ts);
Vfbarra = pontoplanoz(zeta,double(xwn(zeta,Ts)),Ts);

pontos1 = [0 pi/(1.5*sqrt(1-zeta^2)*Ts)];
pontos2 = pontos1;
vec1 = [Vo Vfbarra];
vec2 = vec1;

Satual = polyshape(real([vec2 Vf]),imag([vec2 Vf])).area;

Fp = [Fp, (A*P+B*Z+Z'*B'+P*A'-2*real(Vf)*P>=0):'Limitação à direita'];
Fp = [Fp, (lmisetorconico(loc(Vo,Vfbarra),acos(abs(Vo-real(Vfbarra))/abs(Vo-Vfbarra)), ...
  sysd,P,Z,'E')<=0):'Setor cônico Esquerdo'];

sol = optimize(Fp,trace(P),options);

[primalres,dualres] = check(Fp);
primalres = sort(primalres,'ascend');
dualres = sort(dualres,'ascend');

hold on
plot(real(pontoplanoz(zeta,wnv,Ts)),imag(pontoplanoz(zeta,wnv,Ts)));
pgon = polyshape([real([real(Vf) vec2])],[imag([0 vec2])]);
plot(pgon);

while primalres(1) < 0 || dualres(1) < 0
  if k < length(pontos1)-1
    k = k+1;
  else
    k = 1;
    vec1 = vec2;
    pontos1 = pontos2;
  end
  
  Fp = [];
  Fp = [Fp, (P>=0):'Positividade'];
  Fp = [Fp, (lmiestabilidade(sigma,sysd,Ts,P,Z)<=0):'Estabilidade relativa'];

  Vnew = pontoplanoz(zeta,(pontos1(k)+pontos1(k+1))/2,Ts);
  pontos2 = sort([pontos2 (pontos1(k)+pontos1(k+1))/2],'descend');
  vec2 = sort([vec2 Vnew],'descend');

  Fp = [Fp, (A*P+B*Z+Z'*B'+P*A'-2*real(Vf)*P>=0):'Limitação à direita'];
  
  for m=1:length(vec1)-1
    vLoc = loc(vec2(m),vec2(m+1));
    if vLoc < 0
      phi = acos(abs(real(vec2(m+1))-vLoc)/abs(vec2(m+1)-vLoc));
      Fp = [Fp, (lmisetorconico(vLoc,phi,sysd,P,Z,'D')<=0):['Setor cônico direito ' num2str(k)]];
    else
      phi = acos(abs(vLoc-real(vec2(m+1)))/abs(vLoc-vec2(m+1)));
      Fp = [Fp, (lmisetorconico(vLoc,phi,sysd,P,Z,'E')<=0):['Setor cônico esquerdo ' num2str(k)]];
    end
  end

  sol = optimize(Fp,trace(P),options);
  [primalres,dualres] = check(Fp);
  primalres = sort(primalres,'ascend');
  dualres = sort(dualres,'ascend');
  
  Sant = Satual;
  Satual = polyshape(real(vec2),imag(vec2)).area;

  if (Sant/Satual) < 1 && (Sant/Satual > 0.999999)
    disp('Infactível!');
    break;
  end
  pgon = polyshape([real([real(Vf) vec2])],[imag([0 vec2])]);
  plot(pgon);
end

hold on
axis equal
xline(real(Vf),'m')
plot(real(taxadedecaimento(sigma,Ts)),imag(taxadedecaimento(sigma,Ts)),'m')
K = value(Z)/value(P);
syscomp = ss(A+B*K,B,C+D*K,D,Ts);
pzmap(syscomp,'r')
zgrid(zeta,-1,-1)