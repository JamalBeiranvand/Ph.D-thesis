clear
clc
Eu=10^(20/10);
Rh=100;
v=3.8;
vr=(10^(8/10))^2;
meen=2.2;
K=10;
n=10^3;
m=20:40:500;
mu=log(meen^2/sqrt(vr+meen^2));
sig=sqrt(log(vr/(meen^2)+1));
z=lognrnd(mu,sig,K,n);
Rk=unifrnd(100,925,K,n);
 B=z./((Rk./100).^v);
 
 %.............perfect mrc
for l=1:length(m);
    p=Eu/m(l);
for i=1:K
    sum_B=sum(B);
    num=p*(m(l)-1)*B(i,:);
    den=p*(sum_B-B(i,:))+1;
    x=num./den;
    c(i,:)=log2(1+x);
end 
Rp(l)=sum(sum(c)/n);
end 
plot(m,Rp,'--sr')
axis([min(m) max(m) 0 50])
clear sum_B num den c Rp
%..................ZF
for l=1:length(m);
     p=Eu/m(l);
    x=p*(m(l)-K)*B;
    c=log2(1+x);
    Rp(l)=sum(sum(c))/n;
end 
Rp(l)=sum(sum(c)/n);
hold on
plot(m,Rp,'--*')
axis([min(m) max(m) 0 50]);
clear x c Rp
%................mmse
 for l=1:length(m)
      p=Eu/m(l);
 M=m(l);
x=10*ones(K,n);
a=10*ones(K,n);
t=0;
while sum(sum(a))>1
num=p.*B.*x+(1/(K-1));
den=(M*p*B.*(1-((K-1)*(x-1)./M))+1).^2;
sam=num./den;
for i=1:K
x1(i,:)=sum(sam)-sam(i);
end 
 a=abs(x-x1);
 x=x1;
 t=t+1;
end 
muu=x;
 s=(M.*p.*(1-((K-1).*(muu-1)./M)).*B)+1;
 for i=1:K
 y=1./s(i,:);
num1(i,:)=(1/(K-1))*(sum(1./s)-y);
num2(i,:)=1+p*(sum(B./s.^2)-B(i,:).*y.^2);
kisi(i,:)=muu(i,:)./(num1(i,:).*num2(i,:));
 end
alpha=(M-K+1+(K-1)*muu).^2./(M-K+1+(K-1).*kisi);
teta=p*B.*(M-K+1+(K-1).*kisi)./(M-K+1+(K-1)*muu);
clear x
x=(alpha-1).*teta;
c=log2(1+x);
Rp(l)=sum(sum(c)/n);
 end 
 hold on
 plot(m,Rp,'--g')




