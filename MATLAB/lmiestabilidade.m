<<<<<<< HEAD
function LMI = lmiestabilidade(sigma,sys,Ts,P,Z)
  a11 = -exp(-abs(sigma)*Ts)*P;
  a21 = (P*sys.A'+Z'*sys.B');
  
  LMI = [a11 a21';
    a21 a11];
=======
function LMI = lmiestabilidade(sigma,sys,Ts,P,Z)
  a11 = -exp(-abs(sigma)*Ts)*P;
  a21 = (P*sys.A'+Z'*sys.B');
  
  LMI = [a11 a21';
    a21 a11];
>>>>>>> 1a0eb9ee5d9d76c922a4e1fca1f47741208c87a1
end