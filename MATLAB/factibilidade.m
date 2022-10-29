function K = factibilidade(SYS,TS,SIGMA,ZETA,WN,METODO,PLOTAR)
  % FACTIBILIDADE determina se há uma matriz de ganho K capaz de estabilizar o
  % sistema informado com os parâmetros desejados.
  %
  % K = FACTIBILIDADE(SYS,TS,SIGMA,ZETA,WN,METODO,PLOTAR) verifica se um sistema
  % discreto representado via espaço de estados é factível dentro dos parâmetros
  % de projeto informados. São eles:
  %
  %   - SYS       representação do modelo via espaço de estados;
  %   - TS        período de amostragem;
  %   - SIGMA     valor de estabilidade relativa;
  %   - ZETA      taxa de amortecimento;
  %   - WN        frequência natural-amortecida;
  %   - METODO    (OPCIONAL) método de aproximação usada, onde:
  %       - 'C'     aproximação cônica (VALOR-PADRÃO)
  %       - 'E'     aproximação elíptica
  %       - 'P'     aproximação poligonal
  %   - PLOTAR    (OPCIONAL) opção booleana que plota as regiões
  %               aproximadas. Valor-padrão é FALSO.
  arguments
    SYS
    TS
    SIGMA
    ZETA
    WN
    METODO {mustBeMember(METODO,['C','E','P'])} = 'C'
    PLOTAR {mustBeNumericOrLogical} = false
  end
  A = SYS.A;
  B = SYS.B;
  C = SYS.C;
  D = SYS.D;
  [n,m] = size(B);
  
  NY = 2*pi/(WN*TS);
  
  if NY < 4.86
    error('Curva NY não convexa.');
  end

  P = sdpvar(n,n);
  Z = sdpvar(m,n,'f');
  %   mu = sdpvar(1);
  options = sdpsettings('verbose',0);

  F = [];
  F = [F, (P>=0):'Positividade']; %#ok<*BDSCA>

  %% Estabilidade Relativa
  F = [F, (lmiestabilidade(SIGMA,SYS,TS,P,Z)<=0):'Taxa de amortecimento'];
  
  %% Verifica o método escolhido para a aproximação da taxa de amortecimento
  switch METODO
    case 'C'   % Aproximação Cônica da Taxa de Amortecimento
      Vo = pontoplanoz(ZETA,0,TS);
      Vi = pontoplanoz(ZETA,pi/(sqrt(1-ZETA^2)*TS),TS);
      V = determinarmaiorarea(Vo,Vi,ZETA,TS);
      theta1 = acos(abs(real(V)-Vo)/abs(V-Vo));
      theta2 = acos(abs(real(V)-Vi)/abs(V-Vi));
      
      F = [F, (lmisetorconico(real(Vo),theta1,SYS,P,Z,'E')<=0):['Setor cônico esquerdo ZETA']];
      F = [F, (lmisetorconico(real(Vi),theta2,SYS,P,Z,'D')<=0):['Setor cônico direito ZETA']];
      F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Vi)*P>=0):['Limitação à direita ZETA']];

      No = pontoplanoz(0,WN,TS);
      Ni = pontoplanoz(1,WN,TS);
      theta = acos(abs(real(No)-Ni)/abs(No-Ni));
      
      F = [F, (lmisetorconico(Ni,theta,SYS,P,Z,'D')<=0):['Setor cônico direito WN']];
      F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Ni)*P>=0):['Limitação à direita WN']];
    
    case 'E'  % Aproximação Elíptica
      Vo = pontoplanoz(ZETA,0,TS);
      Vi = pontoplanoz(ZETA,pi/(sqrt(1-ZETA^2)*TS),TS);
      V = determinarmaiorarea(Vo,Vi,ZETA,TS);
      theta1 = acos(abs(real(V)-Vo)/abs(V-Vo));
      theta2 = acos(abs(real(V)-Vi)/abs(V-Vi));
      
      F = [F, (lmisetorconico(real(Vo),theta1,SYS,P,Z,'E')<=0):['Setor cônico esquerdo ZETA']];
      F = [F, (lmisetorconico(real(Vi),theta2,SYS,P,Z,'D')<=0):['Setor cônico direito ZETA']];
      F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Vi)*P>=0):['Limitação à direita ZETA']];

      a = 1-exp((-2*pi)/NY);
      b = (a^2*sin((-2*pi)/NY))/(sqrt(a^2-(cos((-2*pi)/NY)-1)^2));
      
      e11 = -P;
      e21 = -1/a*P+1/2*(1/a+1/b)*(A*P+B*Z)+1/2*(1/a-1/b)*(P*A'+Z'*B');
      
      E = [e11 e21';
        e21 e11];

      F = [F, (E<=0):'Elipse'];

    case 'P'  % Aproximação Poligonal
      l = 0;
      Vo = pontoplanoz(ZETA,0,TS);
      Vi = pontoplanoz(ZETA,pi/(sqrt(1-ZETA^2)*TS),TS);
      V = double(pontoplanoz(ZETA,realwn(ZETA,TS),TS));
      No = pontoplanoz(0,WN,TS);
      Ni = pontoplanoz(1,WN,TS);
      
      pontos1 = [0 pi/(1.5*sqrt(1-ZETA^2)*TS)];
      pontos2 = pontos1;
      vec1 = [Vo V];
      vec2 = vec1;
      pontos3 = [0 1];
      pontos4 = pontos3;
      vec3 = [Ni No];
      vec4 = vec3;
      theta = acos(abs(real(vec4(m+1))-Ni)/abs(vec4(m+1)-Ni));
      
      Satual = polyshape(real([vec2 Vi]),imag([vec2 Vi])).area;
      
      F = [];
      F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Vi)*P>=0):['Limitação à direita ZETA ' METODO]];
      F = [F, (lmisetorconico(loc(Vo,V),acos(abs(Vo-real(V))/abs(Vo-V)), ...
        SYS,P,Z,'E')<=0):['Setor cônico esquerdo ZETA' METODO]];
      F = [F, (lmisetorconico(Ni,theta,SYS,P,Z,'D')<=0):['Setor cônico direito NY ' METODO]];
      F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Ni)*P>=0):['Limitação à direita NY ' METODO]];
      
      optimize(F,trace(P),options);
      
      [primalres,dualres] = check(F);
      primalres = sort(primalres,'ascend');
      dualres = sort(dualres,'ascend');
      
      while primalres(1) < 0 || dualres(1) < 0
        if l < length(vec1)-1
          l = l+1;
        else
          l = 1;
          vec1 = vec2;
          pontos1 = pontos2;
          vec3 = vec4;
          pontos3 = pontos4;
        end
        
        F = [];
        F = [F, (P>=0):'Positividade'];
        F = [F, (lmiestabilidade(SIGMA,SYS,TS,P,Z)<=0):['Taxa de amortecimento ' METODO]];

        Vnew1 = pontoplanoz(ZETA,(pontos1(l)+pontos1(l+1))/2,TS);
        pontos2 = sort([pontos2 (pontos1(l)+pontos1(l+1))/2],'descend');
        vec2 = sort([vec2 Vnew1],'descend');

        Vnew2 = pontoplanoz((pontos3(l)+pontos3(l+1))/2,WN,TS);
        pontos4 = sort([pontos4 (pontos3(l)+pontos3(l+1))/2],'ascend');
        vec4 = sort([vec4 Vnew2],'descend');
      
        F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Vi)*P>=0):['Limitação à direita ZETA ' METODO]];
        F = [F, (A*P+B*Z+Z'*B'+P*A'-2*real(Ni)*P>=0):['Limitação à direita NY ' METODO]];
        
        for m=1:length(vec1)-1
          vLoc = loc(vec2(m),vec2(m+1));
          if vLoc < 0
            phi = acos(abs(real(vec2(m+1))-vLoc)/abs(vec2(m+1)-vLoc));
            F = [F, (lmisetorconico(vLoc,phi,SYS,P,Z,'D')<=0):['Setor cônico direito ZETA ' METODO ' ' num2str(l)]];
          else
            phi = acos(abs(vLoc-real(vec2(m+1)))/abs(vLoc-vec2(m+1)));
            F = [F, (lmisetorconico(vLoc,phi,SYS,P,Z,'E')<=0):['Setor cônico esquerdo ZETA ' METODO ' ' num2str(l)]];
          end
        end

        for m=1:length(vec3)-1
          vLoc2 = loc(vec4(m),vec4(m+1));
          theta = acos(abs(real(vec4(m))-vLoc2)/abs(vec4(m)-vLoc2));
          F = [F, (lmisetorconico(vLoc2,theta,SYS,P,Z,'D')<=0):['Setor cônico direito NY ' METODO ' ' num2str(m)]];
        end
      
        optimize(F,[],options);
        [primalres,dualres] = check(F);
        primalres = sort(primalres,'ascend');
        dualres = sort(dualres,'ascend');
        
        Sant = Satual;
        Satual = polyshape(real(vec2),imag(vec2)).area;
      
        if (Sant/Satual) < 1 && (Sant/Satual > 0.999999)
          break;
        end
      end
  end

  %%
  
  optimize(F,[],options);
  [primalres,dualres] = check(F);
  
  if any(primalres < 0) || any(dualres < 0)
    check(F);
    disp('Infactível!');
  else
    disp('Factível!');
  end

  K = value(Z)/value(P);

  if PLOTAR == true
    epsilon = 1e-3;
    wnv = 0:epsilon:pi/(TS*sqrt(1-ZETA^2));
    zetav = 0:epsilon:1;
    
    cdr = pontoplanoz(ZETA,wnv,TS);
    nfc = pontoplanoz(zetav,WN,TS);
    drc = taxadedecaimento(SIGMA,TS);
  
    hold on
    plot(real(cdr),imag(cdr),':k', ...
      real(cdr),-imag(cdr),':k', ...
      real(nfc),imag(nfc),':k', ...
      real(nfc),-imag(nfc),':k', ...
      real(drc),imag(drc),'--m')
  
    switch METODO
      case 'C'
        pgon = polyshape([real(Vo) real(V) real(Vi) real(V)], ...
          [imag(Vo) imag(V) imag(Vi) -imag(V)]);
        plot(pgon,'LineStyle','--', ...
          'FaceAlpha',0, ...
          'EdgeColor','m')
  
        plot([real(pontoplanoz(0,WN,TS)), exp(-2*pi/NY), real(pontoplanoz(0,WN,TS))], ...
          [imag(pontoplanoz(0,WN,TS)), 0, -imag(pontoplanoz(0,WN,TS))], ...
          'LineStyle','--', ...
          'Color','m')
  
      case 'E'
        pgon = polyshape([real(Vo) real(V) real(Vi) real(V)], ...
          [imag(Vo) imag(V) imag(Vi) -imag(V)]);
        plot(pgon,'LineStyle','--', ...
          'FaceAlpha',0, ...
          'EdgeColor','m')
  
        syms u v
        fimplicit((u-1)^2/a^2+(a^2-(cos((-2*pi)/NY)-1)^2)*v^2/(a^2*sin((2*pi)/NY)^2)==1, ...
          [exp(-2*pi/NY) real(pontoplanoz(0,WN,TS)) -imag(pontoplanoz(0,WN,TS)) imag(pontoplanoz(0,WN,TS))], ...
          'LineStyle','--', ...
          'Color','m')
%         T = [];
%         T(:,1) = fimp.XData';
%         T(:,2) = fimp.YData';
%         writematrix(T, 'data.txt','Delimiter',' ');
 
      case 'P'
        plot(real([vec2 vec2(length(vec2))]),imag([vec2 0]),'--m', ...
          real([vec2 vec2(length(vec2))]),-imag([vec2 0]),'--m')
        plot(real(vec4),imag(vec4),'--m', ...
          real(vec4),-imag(vec4),'--m')
    end
    
    xlabel('Re')
    ylabel('Im')
    SYSCOMP = ss(A+B*K,B,C+D*K,D,TS);
    pzmap(SYSCOMP,'r')
    zgrid(ZETA,WN,TS)
    hold off
  end
end