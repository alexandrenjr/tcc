function LMI = lmiestabilidade(sigma,sys,Ts,P,Z)
  a11 = -exp(-abs(sigma)*Ts)*P;
  a21 = (P*sys.A'+Z'*sys.B');
  
  LMI = [a11 a21';
    a21 a11];