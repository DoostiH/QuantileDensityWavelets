function q=Jons(z,h)
n=length(z);
Q=sort(z);
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
q=1./f1;

