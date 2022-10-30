function V = determinarmaiorarea(Xa,Xb,zeta,Ts)
  fa = 0;
  fb = pi/(Ts*sqrt(1-zeta^2));
  fc = fb/2;
  
  while 1
    L = pontoplanoz(zeta,(fc+fa)/2,Ts);
    N = pontoplanoz(zeta,(fb+fc)/2,Ts);
    
    area1 = 1/2*abs(det([real(Xa),imag(Xa),1; ...
      real(Xb),imag(Xb),1; ...
      real(L),imag(L),1]));
    area2 = 1/2*abs(det([real(Xa),imag(Xa),1; ...
      real(Xb),imag(Xb),1; ...
      real(N),imag(N),1]));
  
    if area1 > area2
      fb = fc;
      fc = (fb+fa)/2;
    elseif area2 > area1
      fa = fc;
      fc = (fb+fa)/2;
    else
      V = L;
      break
    end
  end
end