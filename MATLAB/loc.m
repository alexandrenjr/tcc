function loc = loc(v1,v2)
  if v2 > v1
    loc = real(v1)-imag(v1)/((imag(v2)-imag(v1))/(real(v2)-real(v1)));
  else
    loc = real(v2)-imag(v2)/((imag(v1)-imag(v2))/(real(v1)-real(v2)));
  end
end