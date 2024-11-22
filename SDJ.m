function q=SDJ(z,h)
n=length(z);
u=linspace(1/(n+1),n/(n+1),n);
t=linspace(1/n,1,n);
Q=sort(z);
k=zeros(n,n);
kj=zeros(n,n);
for i=1:n
      kj(:,i)=(1-abs((Q(i)-z)./h));
      for j=1:n
          if kj(j,i)<0
          kj(j,i)=0;
          end
      end
end
f1=sum(kj)/(n*h);
for i=1:n
      k(:,i)=(1-abs((u(i)-t)./h))./f1;
      for j=1:n
          if k(j,i)<0
          k(j,i)=0;
          end
      end
end
q=sum(k)/(n*h);
