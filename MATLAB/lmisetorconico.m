function LMI = lmisetorconico(a,phi,sys,P,Z,direcao)
  arguments
    a
    phi
    sys
    P
    Z
    direcao {mustBeMember(direcao,['E','D'])} = 'E'
  end

  switch direcao
    case 'D'
      a11 = sin(phi)*(2*a*P-sys.A*P-sys.B*Z-P*sys.A'-Z'*sys.B');
      a21 = cos(phi)*(P*sys.A'+Z'*sys.B'-sys.A*P-sys.B*Z);
      
      LMI = [a11 a21';
        a21 a11];

    case 'E'
      a11 = sin(phi)*(sys.A*P+sys.B*Z+P*sys.A'+Z'*sys.B'-2*a*P);
      a21 = cos(phi)*(P*sys.A'+Z'*sys.B'-sys.A*P-sys.B*Z);
      LMI = [a11 a21';
        a21 a11];
  end
end