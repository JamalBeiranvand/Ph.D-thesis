clear
clc
Eu=10^(20/10);
Rh=100;
v=3.8;
vr=(10^(8/10))^2;
meen=2.4;
K=10;
t=K;
T=196;
n=10^3;
m=20:40:500;
mu=log(meen^2/sqrt(vr+meen^2));
sig=sqrt(log(vr/(meen^2)+1));
z=lognrnd(mu,sig,K,n);
Rk=unifrnd(100,925,K,n);
 B=z./((Rk./100).^v);
 
 %.............imperfect mrc
for l=1:length(m);
    p=Eu/(m(l));
for i=1:K
    sum_B=sum(B);
    num=t*p*(m(l)-1)*B(i,:).^2;
    den=(t*p.*B(i,:)+1).*(sum_B-B(i,:))+(t+1).*B(i,:)+1/p;
    x=num./den;
    c(i,:)=log2(1+x);
end 
Rp(l)=((T-t)/T)*sum(sum(c)/n);
end 
plot(m,Rp,'-sr')
axis([min(m) max(m) 0 50])
clear sum_B num den c Rp
 %.............imperfect zf
for l=1:length(m);
    p=Eu/(m(l));
for i=1:K
    sum_B=sum(B./(t*p*B+1));
    num=t*p*(m(l)-K)*B(i,:).^2;
    den=(t*p.*B(i,:)+1).*(sum_B)+(t).*B(i,:)+1/p;
    x=num./den;
    c(i,:)=log2(1+x);
end 
Rp(l)=((T-t)/T)*sum(sum(c)/n);
end 
hold on
plot(m,Rp,'-*')
axis([min(m) max(m) 0 50])
clear sum_B num den c Rp
 %................. imperfect  mmse

 for l=1:length(m)
     p=Eu/(m(l));
      w=1./(sum((B./(t*p*B+1)))+1/p);
 w=repmat(w,K,1);
 B1=(t*p*(B.^2))./(t*p*B+1);
 M=m(l);
x=10*ones(K,n);
a=10*ones(K,n);
while sum(sum(a))>1
num=w.*B1.*x+(1/(K-1));
den=(M*w.*B1.*(1-((K-1)*(x-1)./M))+1).^2;
sam=num./den;
for i=1:K
x1(i,:)=sum(sam)-sam(i);
end 
 a=abs(x-x1);
 x=x1;
end 
muu=x;
 s=(M.*w.*(1-((K-1).*(muu-1)./M)).*B1)+1;
 for i=1:K
 y=1./s(i,:);
num1(i,:)=(1/(K-1))*(sum(1./s)-y);
num2(i,:)=1+w(1,:).*(sum(B1./s.^2)-B1(i,:).*y.^2);
kisi(i,:)=muu(i,:)./(num1(i,:).*num2(i,:));
 end
alpha=(M-K+1+(K-1)*muu).^2./(M-K+1+(K-1).*kisi);
teta=w.*B1.*(M-K+1+(K-1).*kisi)./(M-K+1+(K-1)*muu);
clear x
x=(alpha-1).*teta;
c=log2(1+x);
Rp(l)=((T-t)/T)*sum(sum(c)/n);
 end 
 hold on
 plot(m,Rp,'-g')