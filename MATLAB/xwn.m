function wn = xwn(zeta,Ts)
  syms x
  wn = vpasolve(exp(-Ts*x*cos(acos(zeta)))*cos(Ts*x*sin(acos(zeta)))==exp(-pi*zeta/sqrt(1-zeta^2)));
end