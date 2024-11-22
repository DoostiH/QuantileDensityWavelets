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
N=500;% number of replications
rep1=zeros(n,N);
rep2=zeros(n,N);
rep2=zeros(n,N);

wf = [0.0386 -0.1270 -0.0772 0.6075 0.7457 0.2266]; % For example, filter
% Coiflet(1). 
%parameters in GLD
beta=0.5;
alfa=0.5;
for rep=1:N
%Generating random sample from GLD(landa1,landa2,landa3,landa4)
y=betarnd(alfa,beta,n,1);
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

h = 0.250; % bandwidth
smooth_ll=[];


for j=1:n
  smooth_ll = [smooth_ll loc_lin(j/(n+1), a, linear_est1, h)];
  

end

a1=linspace(1/(n+1),n/(n+1),n);
%f1=(0.5)./sqrt(1-a1);
f1=1./betapdf(betainv(a1,alfa,beta),alfa,beta);
rep1(:,rep)=(f1-smooth_ll).^2;
rep2(:,rep)=(f1-Jons(y,h)).^2;
rep3(:,rep)=(f1-SDJ(y,h)).^2;

rep
end
m1=zeros(1,8);
sd1=zeros(1,8);
m1(1)= mean(rep1(2,:));
s1(1)= std(rep1(2,:));
m1(2)= mean(rep1(21,:));
s1(2)= std(rep1(21,:));
m1(3)= mean(rep1(41,:));
s1(3)= std(rep1(41,:));
m1(4)= mean(rep1(81,:));
s1(4)= std(rep1(81,:));
m1(5)= mean(rep1(121,:));
s1(5)= std(rep1(121,:));
m1(6)= mean(rep1(161,:));
s1(6)= std(rep1(161,:));
m1(7)= mean(rep1(181,:));
s1(7)= std(rep1(181,:));
m1(8)= mean(rep1(199,:));
s1(8)= std(rep1(199,:));

m2=zeros(1,8);
sd2=zeros(1,8);
m2(1)= mean(rep2(2,:));
s2(1)= std(rep2(2,:));
m2(2)= mean(rep2(21,:));
s2(2)= std(rep2(21,:));
m2(3)= mean(rep2(41,:));
s2(3)= std(rep2(41,:));
m2(4)= mean(rep2(81,:));
s2(4)= std(rep2(81,:));
m2(5)= mean(rep2(121,:));
s2(5)= std(rep2(121,:));
m2(6)= mean(rep2(161,:));
s2(6)= std(rep2(161,:));
m2(7)= mean(rep2(181,:));
s2(7)= std(rep2(181,:));
m2(8)= mean(rep2(199,:));
s2(8)= std(rep2(199,:));

m3=zeros(1,8);
s3=zeros(1,8);
m3(1)= mean(rep3(2,:));
s3(1)= std(rep3(2,:));
m3(2)= mean(rep3(21,:));
s3(2)= std(rep3(21,:));
m3(3)= mean(rep3(41,:));
s3(3)= std(rep3(41,:));
m3(4)= mean(rep3(81,:));
s3(4)= std(rep3(81,:));
m3(5)= mean(rep3(121,:));
s3(5)= std(rep3(121,:));
m3(6)= mean(rep3(161,:));
s3(6)= std(rep3(161,:));
m3(7)= mean(rep3(181,:));
s3(7)= std(rep3(181,:));
m3(8)= mean(rep3(199,:));
s3(8)= std(rep3(199,:));
save('d:\Chesneau-Dewan-Doosti\Codes-For-Quantile-Density-Estimation\Table4-3.mat');

    
