<<<<<<< HEAD
function t = taxadedecaimento(sigma,Ts)
  epsilon = 0:0.01:2*pi/Ts;
  t = exp(-abs(sigma)*Ts)*exp(1i.*epsilon.*Ts);
=======
function t = taxadedecaimento(sigma,Ts)
  epsilon = 0:0.01:2*pi/Ts;
  t = exp(-abs(sigma)*Ts)*exp(1i.*epsilon.*Ts);
>>>>>>> 1a0eb9ee5d9d76c922a4e1fca1f47741208c87a1
end