function resultado = realwn(zeta,Ts)
  fun = @(h) cos(h) == -exp(zeta*(h-pi)/(sqrt(1-zeta^2)));
  h0 = pi/(1.72742*Ts*sqrt(1-zeta^2));
  resultado = fzero(fun,h0);
end