clear
clc
k=10;
n=10^3;
P=10;
t=k;
T=196;
Rh=100;
v=3.8;
vr=(10^(8/10))^2;
meen=2.3;
mu=log(meen^2/sqrt(vr+meen^2));
sig=sqrt(log(vr/(meen^2)+1));
m=20:40:500;
%..................mrC 
for l=1:length(m)
    M=m(l);
for N=1:n
    z=lognrnd(mu,sig,k,1);
Rk=unifrnd(100,925,k,1);
 B=z./((Rk./100).^v);
 D=diag(B);
H=sqrt(.5)*(randn(M,k)+1i*randn(M,k));
G=H*sqrt(D);
 A=G;
  A_G=G'*A;
  absA_G=abs((A_G).*conj(A_G));
  num=P*(diag(absA_G)).';
  den=(P*sum(absA_G)-num)+sum(A.*conj(A))*(1+P*sum(B./(t*P*B+1)));
   x(N,:)=num./den;
end 
c=sum(log2(1+x))/n;
Rp(l)=((T-t)/T)*sum(c);
end 
plot(m,Rp,'*')
axis([min(m) max(m) 0 50])
%...............................zf
for l=1:length(m)
    M=m(l);
for N=1:n
    z=lognrnd(mu,sig,k,1);
    Rk=unifrnd(100,925,k,1);
    B=z./((Rk./100).^v);
     D=diag(B);
H=sqrt(.5)*(randn(M,k)+1i*randn(M,k));
G=H*sqrt(D);
 A=G*inv(G'*G);
  A_G=G'*A;
  absA_G=abs((A_G).*conj(A_G));
  %
  num=P*(diag(absA_G)).';
  den=(P*sum(absA_G)-num)+sum(A.*conj(A))*(1+P*sum(B./(t*P*B+1)));
   x(N,:)=num./den;
end 
c=sum(log2(1+x))/n;
Rp(l)=((T-t)/T)*sum(c);
end 
hold on
plot(m,Rp,'r*')
axis([min(m) max(m) 0 50])
%.......................
%........................mmse
for l=1:length(m)
    M=m(l);
for N=1:n
    z=lognrnd(mu,sig,k,1);
Rk=unifrnd(100,925,k,1);
 B=z./((Rk./100).^v);
 D=diag(B);
H=sqrt(.5)*(randn(M,k)+1i*randn(M,k));
G=H*sqrt(D);
 A=G*inv(G'*G+eye(k)./P);
  A_G=G'*A;
  absA_G=abs((A_G).*conj(A_G));
  num=P*(diag(absA_G)).';
  den=(P*sum(absA_G)-num)+sum(A.*conj(A))*(1+P*sum(B./(t*P*B+1)));
   x(N,:)=num./den;
end 
c=sum(log2(1+x))/n;
Rp(l)=((T-t)/T)*sum(c);
end 
hold on
plot(m,Rp,'g*')
axis([min(m) max(m) 0 50])
