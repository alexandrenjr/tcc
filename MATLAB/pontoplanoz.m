function z = pontoplanoz(zeta,wn,Ts)
  % PONTOPLANOZ calcula pontos das curvas da taxa de amortecimento
  % (zeta) e da frequência natural no plano z.
  % 
  % PONTOPLANOZ(ZETA,WN,TS) recebe como parâmetros zeta, a frequência
  % natural não-amortecida e o período de amostragem Ts.
  % 
  % PONTOPLANOZ retorna o ponto z calculado por (6).
  %
  % Dado um ponto no plano s (contínuo) representado por:
  %     s = x + jy                                          (1)
  % a representação no plano z é dada pela transformação
  %     z = exp(sTs)                                        (2)
  % onde Ts é a taxa de amostragem. Substituindo (1) em (2),
  % obtém-se a seguinte relação:
  %     z = exp(xTs)*exp(j*y*Ts)                            (3)
  % Ainda, parte real pode ser representada como
  %     Re(s) = -zeta*wn                                    (4)
  % e a parte imaginária como:
  %     Im(s) = wn*sqrt(1-zeta^2)                           (5)
  % assim, (3) pode ser reescrito como:
  %     z = exp(-zeta*wn*Ts).*exp(j*wn*sqrt(1-zeta^2)*Ts)   (6)
  z = exp(-zeta.*wn.*Ts).*exp(1i.*wn.*sqrt(1-zeta.^2)*Ts);
end