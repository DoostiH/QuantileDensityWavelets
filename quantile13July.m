% Example "A note on the nonparametric estimation
% of a quantile density function by wavelet methods"
% by Christophe Chesneau, Isha Dewan and Hassan Doosti. In this example we
% consider a Generalised Lambda density (GLD)

%  The program calls to:
%           1) Phijk.m -- Daubechies-Lagarias Algorithm. 
%           2) dwtr.m and idwtr.m -- forward and inverse wavelet transform.
%           3) loc_lin.m and kernel.m -- local linear smoother with a
%           Gaussian kernel.
%           4) fijk.m -- itegrate of scale function on [delta-0.5,delta+0.5] 
%           5) simp.m-- estimation of wavelets coefficents with Simpson's
%           Rule
%           6) empiricalcdf.m-- estimation of empirical CDF in fixed points
%           7) Jons.m -- Jones' estimator for quantile density function
%           8) SDJ.m -- Soni-Dewan-Jain's estimator' estimator for quantile density function

close all;
clear all;
lw=2;
set(0, 'DefaultAxesFontSize', 16);
fs=14;
msize=5;
n=200;%sample size
N=1;% number of replications
rep1=zeros(1,N);
rep2=zeros(1,N);
rep2=zeros(1,N);

wf = [0.0386 -0.1270 -0.0772 0.6075 0.7457 0.2266]; % For example, filter
% Coiflet(1). 
%parameters in GLD
landa1=0.5;
landa2=1;
landa3=2;
landa4=6;
for rep=1:N
%Generating random sample from GLD(landa1,landa2,landa3,landa4)
u=rand(n,1);
y=landa1+ (u.^(landa3)-(1-u).^(landa4))/(landa2);

% Apply Daubechies Lagarias to obtain the linear wavelet estimate 



 j0=5; % default coarsest level

 %calculation of \hat{\c}_{(j_0,k_1)}
     ahat=zeros(2^j0,1);
     for k1=0:2^j0-1
                 
             
         ahat(k1+1)=simp(j0,k1,y);
         
     end
     
     
     linear_est1=2^(j0/2)*ahat;  
      
                

a = linspace(0,1, length(linear_est1));

      % How does the linear estimator look like?
         

csw1 = dwtr(linear_est1, 3, wf); % Forward wavelet transformation 

nn = length(csw1);
nn2 = nn/2;
finest = csw1(nn2+1:nn); %finest details
sigma = 1.4826*median(abs(finest-median(finest)));
% sigma1 = 3;
lambda1 = sqrt(  log(n)/n) * sigma;



cswt1 = csw1 .* ( abs(csw1) > 8*lambda1); %hard threshold
smooth_th1 = Idwtr(cswt1, 3, wf); % Inverse wavelet transformation



%--------------------------------------------------------------------------

h = 0.150; % bandwidth
smooth_ll=[];


for j=1:n
  smooth_ll = [smooth_ll loc_lin(j/(n+1), a, linear_est1, h)];
  

end
linspace(1/(n+1),n/(n+1),n);
a1=linspace(1/(n+1),n/(n+1),n);
%f1=(0.5)./sqrt(1-a1);
f1=(landa3*(a1.^(landa3-1))+landa4*((1-a1).^(landa4-1)))/(landa2);

rep1(rep)=norm(f1-smooth_ll);
rep2(rep)=norm(f1-Jons(y,h));
rep3(rep)=norm(f1-SDJ(y,h));

rep
end
ise1=(rep1.^2)*(a1(n)-a1(1))*(a1(2)-a1(1));
ise2=(rep2.^2)*(a1(n)-a1(1))*(a1(2)-a1(1));
ise3=(rep3.^2)*(a1(n)-a1(1))*(a1(2)-a1(1));

a2=linspace(1/(n+1),n/(n+1),n);
%Estimation of MISE for three estimators
mean(ise1)
mean(ise2)
mean(ise3)

%save('d:\Chesneau-Dewan-Doosti\Codes-For-Quantile-Density-Estimation\quantilen500.0.5-1-2-6.mat');
%Estmation of quantile density 

% How does the linear estimators of quantile densities look like?
figure(1);
 plot( a, (landa3*(a.^(landa3-1))+landa4*((1-a).^(landa4-1)))/(landa2),'k-.') % observed f^Y
hold on
plot(a,linear_est1,'b')
title(['Linear wavelet quantile density GLD(' num2str(landa1),',' num2str(landa2),',' num2str(landa3),',' num2str(landa4),'),n=',num2str(n) ])



% How does the linear estimators of quantile densities look like after hard thresholding?
figure(2);
 plot( a, (landa3*(a.^(landa3-1))+landa4*((1-a).^(landa4-1)))/(landa2),'k-.') % observed f^Y
hold on
plot(a, smooth_th1,'r')
title(['Linear wavelet estimator after thresholding-quantile density GLD(' num2str(landa1),',' num2str(landa2),',' num2str(landa3),',' num2str(landa4),'),n=',num2str(n) ])




% How does it look like after local linear smoothing?
figure(3)
 plot( a, (landa3*(a.^(landa3-1))+landa4*((1-a).^(landa4-1)))/(landa2),'k-.') % observed f^Y
hold on
plot((1:n)/(n+1), smooth_ll, 'g')
title(['Smoothed linear wavelet estimator-quantile density GLD(' num2str(landa1),',' num2str(landa2),',' num2str(landa3),',' num2str(landa4),'),n=',num2str(n) ])

% How does it look like three different estimators of quntile density?
    %--------------------------------------------------------------------------
    figure(4)
plot(a, (landa3*(a.^(landa3-1))+landa4*((1-a).^(landa4-1)))/(landa2),'k', a,linear_est1,'b:',a, smooth_th1,'r-.', (1:n)/(n+1), smooth_ll, 'g--',a2, Jons(y,h),'yo',a2,SDJ(y,h),'mx',...
'MarkerSize',2);                
% 'LineWidth',1,...
%                 'MarkerEdgeColor','b',...
%                 'MarkerFaceColor','w',...
%                  'MarkerSize',5);
% plot(a, (landa3*(a.^(landa3-1))+landa4*((1-a).^(landa4-1)))/(landa2),'k',(1:n)/(n+1), smooth_ll, 'g',a2, Jons(y,h),'y',a2,SDJ(y,h),'m',...
%                 'LineWidth',1,...
%                 'MarkerEdgeColor','b',...
%                 'MarkerFaceColor','w',...
%                 'MarkerSize',5);
title(['Estimation of quantile density GLD(' num2str(landa1),',' num2str(landa2),',' num2str(landa3),',' num2str(landa4),'),n=',num2str(n) ])

%print -dpdf fig3.pdf;