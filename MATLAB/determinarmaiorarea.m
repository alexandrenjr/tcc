function M = determinarmaiorarea(L,N,zeta,Ts)
  fa = 0;
  fb = pi/(Ts*sqrt(1-zeta^2));
  fc = fb/2;
  
  while 1
    M1 = pontoplanoz(zeta,(fc+fa)/2,Ts);
    M2 = pontoplanoz(zeta,(fb+fc)/2,Ts);
    
    area1 = 1/2*abs(det([real(L),imag(L),1; ...
      real(N),imag(N),1; ...
      real(M1),imag(M1),1]));
    area2 = 1/2*abs(det([real(L),imag(L),1; ...
      real(N),imag(N),1; ...
      real(M2),imag(M2),1]));
  
    if area1 > area2
      fb = fc;
      fc = (fb+fa)/2;
    elseif area2 > area1
      fa = fc;
      fc = (fb+fa)/2;
    else
      M = M1;
      break
    end
  end
end