function yy = Phijk(z, j, k, filter, n)
%--------------------------------------------------------------
%  yy=Phijk(z, j, k, filter, n)
%  Evaluation of the scaling function corresponding to an Orthogonal 
%           MRA by Daubechies-Lagarias Algorithm.
%  inputs:  z -- the argument
%           j -- scale
%           k -- shift
%           filter -- ON finite wavelet filter, might be an
%                     output of WaveLab's: MakeONFilter
%           n -- precision of approximation maesured by the number
%                of Daubechies-Lagarias steps (default n=20)
%--------------------------------------------------------
%  output:  yy -- value of father wavelet (j,k) coresponding to
%                  'filter' at z.
%--------------------------------------------------------------
% Example of use: 
% > xx = 0:0.01:6;  yy=[]; 
% > for i=1:length(xx)
% > yy =[yy Phijk(x(i), 0, 1, MakeONFilter('Daubechies',4), 25)];
% > end
% > plot(x,yy)
%---------------------------------------------------------------- 
    if (nargin == 4) 
        n=20;
    end
    daun=length(filter)/2;
    N=length(filter)-1;
    x=(2^j)*z-k;
    if(x<=0||x>=N) yy=0;
else
    int=floor(x);
    dec=x-int;
    dy=dec2bin(dec,n);
    t0=t0(filter);
    t1=t1(filter);
    prod=eye(N);
    for i=1:n
            if dy(i)==1 prod=prod*t1;
            else prod=prod*t0;
            end
    end
    y=2^(j/2)*prod;
      yyy = mean(y');
      yy = yyy(int+1);
     end
        
  %--------------------functions needed----------------------------
  %---------
function a = dec2bin(x,n)
            a=[];
             for i = 1:n
             if(x <= 0.5) a=[a 0]; x=2*x;
              else a=[a 1]; x=2*x-1;
             end
            end
 %-----------
 function t0 = t0(filter)
%
n = length(filter);
nn = n - 1;
%
t0 = zeros(nn);
for i = 1:nn
    for j= 1:nn
        if (2*i - j > 0 & 2*i - j <= n) 
            t0(i,j) = sqrt(2) * filter( 2*i - j ); 
        end
    end
end
%------------------
function t1 = t1(filter)
%
n = length(filter);
nn = n - 1;
%
t1 = zeros(nn);
for i = 1:nn
    for j= 1:nn
        if (2*i -j+1 > 0 & 2*i - j+1 <= n) 
            t1(i,j) = sqrt(2) * filter( 2*i - j+1 ); 
        end
    end
end
%---------------------- B. Vidakovic, 2002 --------------------
%-----------------------------------------------------------------