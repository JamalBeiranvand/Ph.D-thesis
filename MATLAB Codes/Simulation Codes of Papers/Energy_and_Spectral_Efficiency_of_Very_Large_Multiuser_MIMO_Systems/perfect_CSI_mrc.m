clear
clc
p=10;
Rh=100;
v=3.8;
vr=(10^(8/10))^2;
meen=2.5;
K=10;
n=1000;
m=20:40:500;
mu=log(meen^2/sqrt(vr+meen^2));
sig=sqrt(log(vr/(meen^2)+1));
z=lognrnd(mu,sig,K,n);
Rk=unifrnd(100,925,K,n);
 B=z./((Rk./100).^v);
 for l=1:length(m)
     M=m(l);
  H=sqrt(.5)*(randn(M*K,n)+1j*randn(M*K,n));   
 for i=1:K
     t=((i-1)*M)+(1:1:M);
    D(t,:)=(repmat(B(i,:),M,1));
 end 
 G=H.*sqrt(D);
 A=conj(G);
 
 for i=1:K
     t=((i-1)*M)+(1:1:M);
     a=A(t,:);
     for j=1:K
       re=((j-1)*M)+(1:1:M);  
       s(j,:)=p.*abs(sum(a.*G(re,:)));
     end 
       num=s(i,:);
       den=sum(s)-num+(sum((a.*conj(a))));
       x=num./den;
       c(i,:)=log2(1+x);
 end 
    Rp(l)=sum(sum(c))/n ;
 end
 hold on
 plot(m,Rp)
axis([min(m) max(m) 0 50])
         