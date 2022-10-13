function t = taxadedecaimento(sigma,Ts)
  epsilon = 0:0.01:2*pi/Ts;
  t = exp(-abs(sigma)*Ts)*exp(1i.*epsilon.*Ts);
end