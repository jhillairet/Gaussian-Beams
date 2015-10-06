function s = csgn(x)
%  The csgn function is used to determine in which half-plane ("left" or
%    "right") the complex-valued expression or number x lies. It is defined by
%  
%                 /  1    if Re(x) > 0 or Re(x) = 0 and Im(x) > 0
%      csgn(x) = <
%                 \ -1    if Re(x) < 0 or Re(x) = 0 and Im(x) < 0
if real(x) ~= 0
  s = sign(real(x));
else
  s = sign(imag(x));
end