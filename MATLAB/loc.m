<<<<<<< HEAD
function loc = loc(v1,v2)
  if v2 > v1
    loc = real(v1)-imag(v1)/((imag(v2)-imag(v1))/(real(v2)-real(v1)));
  else
    loc = real(v2)-imag(v2)/((imag(v1)-imag(v2))/(real(v1)-real(v2)));
  end
=======
function loc = loc(v1,v2)
  if v2 > v1
    loc = real(v1)-imag(v1)/((imag(v2)-imag(v1))/(real(v2)-real(v1)));
  else
    loc = real(v2)-imag(v2)/((imag(v1)-imag(v2))/(real(v1)-real(v2)));
  end
>>>>>>> 1a0eb9ee5d9d76c922a4e1fca1f47741208c87a1
end