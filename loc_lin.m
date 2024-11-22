function smoo=loc_lin(x, vecX, vecY, bandwith)
% Local linear regression
% Usage
%   smoo=loc_lin(x, vecX, vecY, bandwith)
% Input
%   x, vecX, vecY, bandwith
% Output
%   smoo
% See also
%   kernel, nada_wat
% Brani  11/2002

N=length(vecX);
  num   = 0; denom = 0; s1=0; s2=0;
  for i=1:N
  s1 = s1 + kernel(vecX(i) - x, bandwith)*(vecX(i)-x);
  s2 = s2 + kernel(vecX(i) - x, bandwith)*(vecX(i)-x)^2;
  end
%---------
for i = 1  :  N
  wei(i) = kernel(vecX(i) - x, bandwith) * (s2 - (vecX(i) - x)* s1);
  num=num+wei(i)*vecY(i);
  denom=denom+wei(i);
end
%==================
smoo = num/denom;

