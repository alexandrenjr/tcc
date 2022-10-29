clear;
clc;

epsilon = 1e-6;
zeta = 0.9;
Ts = 1;
wnmax = pi/(sqrt(1-zeta^2)*Ts);
wn = wnmax/2;

e = pontoplanoz(zeta,wnmax,Ts);
f = pontoplanoz(zeta,wn,Ts);

while real(e) ~= real(f)
  if real(e) < real(f)
   wn = wn + epsilon;
  else
    wn = wn - epsilon;
  end

  if abs(real(e) - real(f)) < epsilon
    break;
  end

  f = pontoplanoz(zeta,wn,Ts);
end

hold
plot([real(e),real(f)],[imag(e),imag(f)])
zgrid(zeta,-1,-1)