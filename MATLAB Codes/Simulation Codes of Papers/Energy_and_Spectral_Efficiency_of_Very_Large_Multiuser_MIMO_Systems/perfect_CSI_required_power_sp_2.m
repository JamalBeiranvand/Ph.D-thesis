clear
clc
sp=2;
K=10;
Rp=sp*K;
Rh=100;
v=3.8;
vr=(10^(8/10))^2;
meen=2.2;
n=10^3;
mu=log(meen^2/sqrt(vr+meen^2));
sig=sqrt(log(vr/(meen^2)+1));
z=lognrnd(mu,sig,K,n);
Rk=unifrnd(100,925,K,n);
 B=z./((Rk./100).^v);

 %.............perfect mrc
m=100:40:500;
power(1)=26.9
 p=30;
for l=2:length(m);
    a=1
    while a>10^-1
        p=p-.001;
for i=1:K
    sum_B=sum(B);
    num=p*(m(l)-1)*B(i,:);
    den=p*(sum_B-B(i,:))+1;
    x=num./den;
    c(i,:)=log2(1+x);
end 
Rp1=sum(sum(c)/n);
a=abs(Rp-Rp1);
    end 
    power(l)=10*log10(p);
end 
plot(m,power,'*-')
axis([min(m) max(m) -3 30])
clear sum_B num den c
%..................ZF
p=60;
m=20:40:500;
for l=1:length(m);
    a=2
    while a>10^-2
        p=p-.001;
    x=p*(m(l)-K)*B;
    c=log2(1+x);
    Rp1=sum(sum(c))/n;
    a=abs(Rp-Rp1);
    end 
    power(l)=10*log10(p);
    
end 
hold on
plot(m,power,'>-')
axis([min(m) max(m) -3 30])
clear sum_B num den c 
%................mmse
m=20:40:500;
 p=50;
 for l=1:length(m)
 M=m(l);
dif=1
while dif >.1
     x=10*ones(K,n);
a=10*ones(K,n);
    p=p-.005;
    w=0;
while sum(sum(a))>1 & w<100
    w=w+1;
num=p.*B.*x+(1/(K-1));
den=(M*p*B.*(1-((K-1)*(x-1)./M))+1).^2;
sam=num./den;
for i=1:K
x1(i,:)=sum(sam)-sam(i);
end 
 a=abs(x-x1);
 x=x1;
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
Rp1=sum(sum(c)/n);
dif=abs(Rp-Rp1);
    end 
    power(l)=10*log10(p);
 end 
 hold on
plot(m,power,'-')
axis([min(m) max(m) -3 30])
clear sum_B num den c