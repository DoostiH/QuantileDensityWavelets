function y = kernel(x, bandwith)
% Gaussian kernel
% Usage
%   y = kernel(x, bandwith)
% Input
%   x, bandwith
% Output
%   y
% See also
%   nada_wat, loc_lin
% Brani  11/2002

y = 1/bandwith * 1/sqrt(2 * pi) * exp( - 1/2 * (x/bandwith).^2 );
