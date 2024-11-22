function  data = idwtr(wtr, L, filterh)
% function data = idwt(wtr, L, filterh); Calculates the IDWT of wavelet
% transformation wtr using wavelet filter  "filterh"  and  L  scales.  
% Use
%>> max(abs(data - IDWTR(DWTR(data,3,filter), 3,filter)))
%
%ans = 4.4409e-016

nn = length(wtr);   n = length(filterh);           % Lengths
if nargin==2, L = round(log2(nn)); end;            % Depth of transformation
H = filterh;                                       % Wavelet H filter
G = fliplr(H); G(2:2:n) = -G(2:2:n);               % Wavelet G filter
LL = nn/(2^L);                                     % Number of scaling coeffs
C =  wtr(1:LL);                                    % Scaling coeffs
for j = 1:L                                        % Cascade algorithm
   w  = mod(0:n/2-1,LL)+1;                         % Make periodic
   D  = wtr(LL+1:2*LL);                            % Wavelet coeffs
   Cu(1:2:2*LL+n) = [C C(1,w)];                    % Upsample & keep periodic
   Du(1:2:2*LL+n) = [D D(1,w)];                    % Upsample & keep periodic
   C  = conv(Cu,H) + conv(Du,G);                   % Convolve & add
   C  = C([n:n+2*LL-1]-1);                         % Periodic part
   LL = 2*LL;                                      % Double the size of level
end;
data = C;                                          % The inverse DWT

