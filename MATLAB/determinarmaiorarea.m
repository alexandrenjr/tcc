function Q = determinarmaiorarea(Xa,Xb,zeta,Ts)
  fa = 0;
  fb = pi/(Ts*sqrt(1-zeta^2));
  fc = fb/2;
  
  while 1
    Q1 = pontoplanoz(zeta,(fc+fa)/2,Ts);
    Q2 = pontoplanoz(zeta,(fb+fc)/2,Ts);
    
    area1 = 1/2*abs(det([real(Xa),imag(Xa),1; ...
      real(Xb),imag(Xb),1; ...
      real(Q1),imag(Q1),1]));
    area2 = 1/2*abs(det([real(Xa),imag(Xa),1; ...
      real(Xb),imag(Xb),1; ...
      real(Q2),imag(Q2),1]));
  
    if area1 > area2
      fb = fc;
      fc = (fb+fa)/2;
    elseif area2 > area1
      fa = fc;
      fc = (fb+fa)/2;
    else
      Q = Q1;
      break
    end
  end