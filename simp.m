function y = simp (j0,k,sam )

%SIMP       approximate the definite integral of an arbitrary function
%           using the composite Simpson's rule (i.e., the closed
%           Newton-Cotes quadrature rule with n=2)
%
%     calling sequences:
%             y = simp ( 'f', a, b, n )
%             simp ( 'f', a, b, n )
%
%     inputs:
%             f       string containing name of m-file defining integrand
%             a       lower limit of integration
%             b       upper limit of integration
%             n       number of uniformly sized subintervals into which
%                     integration interval is to be divided - must be even
%                     (the resulting approximation will require n+1
%                     function evaluations)
%
%     output:
%             y       approximate value of the definite integral of f(x)
%                     over the interval a < x < b
%
%     NOTE:
%             if SIMP is called with no output arguments, the approximate 
%             value of the definite integral of f(x) over the interval 
%             a < x < b will be displayed
%
%To try another filter: wf = MakeONFilterExt(Type,Par)
% for example Symmlet 4:

filter = [   0.332670552950082615998512,  0.806891509311092576494494, ...
                0.459877502118491570095152, -0.135011020010254588696390, ...
               -0.085441273882026661692819,  0.0352262918857095366027407 	];
% filter=[ 	0.2303778133088965,     0.7148465705529158, ...
%                 0.630880767929859,     -0.02798376941686011, ...
%                -0.1870348117190932,     0.0308413818355608, ...
%                 0.03288301166688522,   -0.01059740178506904 	];
             wf = [0.0386 -0.1270 -0.0772 0.6075 0.7457 0.2266]; % For example, filter
% Coiflet(1). 
n=length(sam);
%b=max(sam)+2*std(sam);
% b=landa1+1/(landa2);
% a=landa1-1/(landa2);
b=max(sam);
a=min(sam);
h = (b-a)/n;
x = linspace ( a, b, n+1 );
t=empiricalcdf(sam,x);
% t(1)=0;
% t(2)=0;
j3=find(t==min(t), 1 );
for i = 1:j3-1
    t(i)=0;
end
j1=find(t==max(t), 1 );
for i = j1:length(t)
    t(i)=1;
end
for i = 1:n+1
    fx(i) = Phijk(t(i),j0,k,filter,25);
end
w = [ 1 zeros(1,n-1) 1 ];
w(2:2:n) = 4*ones(1,n/2);
w(3:2:n-1) = 2*ones(1,n/2-1);

if ( nargout == 1 ) 
   y = (h/3) * sum ( w .* fx );
else
   disp ( (h/3) * sum ( w .* fx ) );
end