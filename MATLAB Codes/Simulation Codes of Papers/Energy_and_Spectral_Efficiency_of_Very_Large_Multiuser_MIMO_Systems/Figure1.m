clear
clc
p=10;
Rh=100;
v=3.8;
vr=(10^(8/10));
K=10;
n=10^3;
MM=20:40:500;
z=exp(sqrt(vr)*sqrt(.5)*(randn(K,n)+1j*randn(K,n)));
Rk=unifrnd(Rh,1000,K,n);
B=z./((Rk./Rh).^v);
 
 %.............perfect mrc
for l=1:length(MM)
    for i=1:K
        sum_B=sum(B);
        num=p*(MM(l)-1)*B(i,:);
        den=p*(sum_B-B(i,:))+1;
        x=num./den;
        c(i,:)=log2(1+x);
    end 
    Rp(l)=sum(sum(c)/n);
end 
subplot(2,1,1)
plot(MM,Rp)
axis([min(MM) max(MM) 0 50])
grid on
clear sum_B num den c Rp
% %..................ZF
for l=1:length(MM)
    x=p*(MM(l)-K)*B;
    c=log2(1+x);
    Rp(l)=sum(sum(c))/n;
end 
Rp(l)=sum(sum(c)/n);
subplot(2,1,1)
hold on
plot(MM,Rp,'r')
axis([min(MM) max(MM) 0 50]);
clear x c Rp
%................mmse
 for l=1:length(MM)
    M=MM(l);
    x=10*ones(K,n);
    a=10*ones(K,n);
    
    while sum(sum(a))>K*n
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
    Rp(l)=sum(sum(c)/n);
  end 
 subplot(2,1,1)
 hold on
 plot(MM,Rp,'g')
% %.........................................
clear
clc
K=10;
n=10^3;
P=10;
Rh=100;
v=3.8;
vr=(10^(8/10));
MM=20:40:500;
%..................mrC 
for l=1:length(MM)
    M=MM(l);
    for N=1:n
        z=exp(sqrt(vr)*sqrt(.5)*(randn(K,1)+1j*randn(K,1)));
        Rk=unifrnd(Rh,1000,K,1);
        B=z./((Rk./Rh).^v);
        D=diag(B);
        H=sqrt(.5)*(randn(M,K)+1i*randn(M,K));
        G=H*sqrt(D);
        A=G;
        A_G=G'*A;
        absA_G=abs((A_G).*conj(A_G));
        num=P*(diag(absA_G)).';
        den=P*sum(absA_G)-num+sum(A.*conj(A));
        x(N,:)=num./den;
    end 
    c=sum(log2(1+x))/n;
    Rp(l)=sum(c);
end 
hold on
subplot(2,1,1)
plot(MM,Rp,'b*')
axis([min(MM) max(MM) 0 50])
%.........................ZF
for l=1:length(MM)
    M=MM(l);
    for N=1:n
            z=exp(sqrt(vr)*sqrt(.5)*(randn(K,1)+1j*randn(K,1)));
            Rk=unifrnd(Rh,1000,K,1);
            B=z./((Rk./100).^v);
            D=diag(B);
            H=sqrt(.5)*(randn(M,K)+1i*randn(M,K));
            G=H*sqrt(D);
            A=G*inv(G'*G);
            A_G=G'*A;
            absA_G=abs((A_G).*conj(A_G));
            num=P*(diag(absA_G)).';
            den=P*sum(absA_G)-num+sum(A.*conj(A));
            x(N,:)=num./den;
    end 
    c=sum(log2(1+x))/n;
    Rp(l)=sum(c);
end 
hold on
subplot(2,1,1)
plot(MM,Rp,'r*')
axis([min(MM) max(MM) 0 50])
%........................mmse
for l=1:length(MM)
    M=MM(l);
    for N=1:n
        z=exp(sqrt(vr)*sqrt(.5)*(randn(K,1)+1j*randn(K,1)));
        Rk=unifrnd(Rh,1000,K,1);
        B=z./((Rk./100).^v);
        D=diag(B);
        H=sqrt(.5)*(randn(M,K)+1i*randn(M,K));
        G=H*sqrt(D);
        A=G*inv(G'*G+eye(K)./P);
        A_G=G'*A;
        absA_G=abs((A_G).*conj(A_G));
        num=P*(diag(absA_G)).';
        den=P*sum(absA_G)-num+sum(A.*conj(A));
        x(N,:)=num./den;
    end 
    c=sum(log2(1+x))/n;
    Rp(l)=sum(c);
end 
hold on
subplot(2,1,1)
plot(MM,Rp,'g*')
axis([min(MM) max(MM) 0 50])
%% ...........................imperfect
clear
clc
p=10;
Rh=100;
v=3.8;
vr=(10^(8/10));
K=10;
t=K;
T=196;
n=10^3;
MM=20:40:500;
z=exp(sqrt(vr)*sqrt(.5)*(randn(K,n)+1j*randn(K,n)));
Rk=unifrnd(100,1000,K,n);
B=z./((Rk./Rh).^v);
 
 %.............imperfect mrc
for l=1:length(MM)
    for i=1:K
        sum_B=sum(B);
        num=t*p*(MM(l)-1)*B(i,:).^2;
        den=(t*p.*B(i,:)+1).*(sum_B-B(i,:))+(t+1).*B(i,:)+1/p;
        x=num./den;
        c(i,:)=log2(1+x);
    end 
    Rp(l)=((T-t)/T)*sum(sum(c)/n);
end 
subplot(2,1,2)
plot(MM,Rp)
clear sum_B num den c Rp
 %...................imperfect zf
for l=1:length(MM)
    for i=1:K
        sum_B=sum(B./(t*p*B+1));
        num=t*p*(MM(l)-K)*B(i,:).^2;
        den=(t*p.*B(i,:)+1).*(sum_B)+(t).*B(i,:)+1/p;
        x=num./den;
        c(i,:)=log2(1+x);
    end 
    Rp(l)=((T-t)/T)*sum(sum(c)/n);
end 
subplot(2,1,2)
hold on
plot(MM,Rp,'r')
clear sum_B num den c Rp
 %.................mmse
 w=1./(sum((B./(t*p*B+1)))+1/p);
 w=repmat(w,K,1);
 B1=(t*p*(B.^2))./(t*p*B+1);
 for l=1:length(MM)
     M=MM(l);
     x=10*ones(K,n);
     a=10*ones(K,n);
     while sum(sum(a))>K*n
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
 subplot(2,1,2)
 hold on
 plot(MM,Rp,'g')
% .................................................
clear
clc
K=10;
n=10^3;
P=10;
t=K;
T=196;
Rh=100;
v=3.8;
vr=(10^(8/10));
MM=20:40:500;
%..................mrC 
for l=1:length(MM)
    M=MM(l);
    for N=1:n
        z=exp(sqrt(vr)*sqrt(.5)*(randn(K,1)+1j*randn(K,1)));
        Rk=unifrnd(Rh,1000,K,1);
        B=z./((Rk./Rh).^v);
        D=diag(B);
        H=sqrt(.5)*(randn(M,K)+1i*randn(M,K));
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
subplot(2,1,2)
plot(MM,Rp,'b*')
%...............................zf
for l=1:length(MM)
    M=MM(l);
    for N=1:n
        z=exp(sqrt(vr)*sqrt(.5)*(randn(K,1)+1j*randn(K,1)));
        Rk=unifrnd(Rh,1000,K,1);
        B=z./((Rk./Rh).^v);
        D=diag(B);
        H=sqrt(.5)*(randn(M,K)+1i*randn(M,K));
        G=H*sqrt(D);
        A=G*inv(G'*G);
        A_G=G'*A;
        absA_G=abs((A_G).*conj(A_G));
        num=P*(diag(absA_G)).';
        den=(P*sum(absA_G)-num)+sum(A.*conj(A))*(1+P*sum(B./(t*P*B+1)));
        x(N,:)=num./den;
    end 
    c=sum(log2(1+x))/n;
    Rp(l)=((T-t)/T)*sum(c);
end 
subplot(2,1,2)
hold on
plot(MM,Rp,'r*')
%.......................
%........................mmse
for l=1:length(MM)
    M=MM(l);
    for N=1:n
        z=exp(sqrt(vr)*sqrt(.5)*(randn(K,1)+1j*randn(K,1)));
        Rk=unifrnd(100,1000,K,1);
        B=z./((Rk./100).^v);
        D=diag(B);
        H=sqrt(.5)*(randn(M,K)+1i*randn(M,K));
        G=H*sqrt(D);
        A=G*inv(G'*G+eye(K)./P);
        A_G=G'*A;
        absA_G=abs((A_G).*conj(A_G));
        num=P*(diag(absA_G)).';
        den=(P*sum(absA_G)-num)+sum(A.*conj(A))*(1+P*sum(B./(t*P*B+1)));
        x(N,:)=num./den;
    end 
    c=sum(log2(1+x))/n;
    Rp(l)=((T-t)/T)*sum(c);
end 
subplot(2,1,2)
hold on
plot(MM,Rp,'g*')
axis([min(MM) max(MM) 0 50])
grid on

