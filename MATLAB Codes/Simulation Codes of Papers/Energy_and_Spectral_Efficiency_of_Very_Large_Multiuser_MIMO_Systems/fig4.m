clear
clc
sp=1;
K=10;
Rp=sp*K;
Rh=100;
v=3.8;
vr=(10^(8/10));
t=K;
T=196;
n=10^3;
z=exp(sqrt(vr)*sqrt(.5)*(randn(K,n)+1j*randn(K,n)));
Rk=unifrnd(Rh,1000,K,n);
 B=z./((Rk./100).^v);
 
 %.............imperfect mrc
 
 m(1)=25;
 p=30;
 power(1)=20.56;
 m(2:13)=60:40:500;
for l=2:length(m)
    
    dif=1;
    while dif>10^-1
          p=p-.01;
          for i=1:K
              sum_B=sum(B);
              num=t*p*(m(l)-1)*B(i,:).^2;
              den=(t*p.*B(i,:)+1).*(sum_B-B(i,:))+(t+1).*B(i,:)+1/p;
              x=num./den;
              c(i,:)=log2(1+x);
          end 
          Rp1=((T-t)/T)*sum(sum(c)/n);
          dif=abs(Rp-Rp1);
     end 
     power(l)=10*log10(p)
     if real(p)<2
         p=10;
     end
end 
1
hold on
plot(m,power,'-*g')
axis([20 max(m) -9 20])
clear sum_B num den c
 %...................imperfect zf
   p=20;
   m=20:40:500;
for l=1:length(m)
    dif=1;
    while dif>10^-1
            p=p-.01;
            for i=1:K
                sum_B=sum(B./(t*p*B+1));
                num=t*p*(m(l)-K)*B(i,:).^2;
                den=(t*p.*B(i,:)+1).*(sum_B)+(t).*B(i,:)+1/p;
                x=num./den;
                c(i,:)=log2(1+x);
            end 
            Rp1=((T-t)/T)*sum(sum(c)/n);
            dif=abs(Rp-Rp1);
    end 
     power(l)=10*log10(p);
     if real(p)<2
         p=5;
     end
end 
2
hold on
plot(m,power,'->g')
axis([min(m) max(m) -9 20])
clear sum_B num den c 
%  %.................mmse
   p=20;
   m=20:40:500;
 for l=1:length(m)
      
        dif=1;
        while dif>10^-1
                p=p-.01;
                w=1./(sum((B./(t*p*B+1)))+1/p);
                w=repmat(w,K,1);
                B1=(t*p*(B.^2))./(t*p*B+1);
                M=m(l);
                x=10*ones(K,n);
                a=10*ones(K,n);
                number=0;
                while sum(sum(a))>K*n && number<100
                      number=number+1;
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
                Rp1=((T-t)/T)*sum(sum(c)/n);
                dif=abs(Rp-Rp1);
        end 
        power(l)=10*log10(p)
 end 
 3
 hold on
 plot(m,power,'-g')
axis([min(m) max(m) -9 20])
%................................perfect
clear
clc
sp=1;
K=10;
Rp=sp*K;
Rh=100;
v=3.8;
vr=(10^(8/10));
n=10^3;
m=20:40:500;
z=exp(sqrt(vr)*sqrt(.5)*(randn(K,n)+1j*randn(K,n)));
Rk=unifrnd(Rh,1000,K,n);
 B=z./((Rk./100).^v);

 %.............perfect mrc
 p=100;
for l=1:length(m)
    a=1;
    while a>10^-1
            p=p-.01;
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
     if real(p)<2
         p=5;
     end
end 
4
hold on
plot(m,power,'*-')
axis([min(m) max(m) -9 20])
clear sum_B num den c

%..................ZF

for l=1:length(m)
    p=12;
    a=2;
    while a>10^-1
        p=p-.001;
        x=p*(m(l)-K)*B;
        c=log2(1+x);
        Rp1=sum(sum(c))/n;
        a=abs(Rp-Rp1);
    end 
    power(l)=10*log10(p);
end 
5
hold on
plot(m,power,'>-')
axis([min(m) max(m) -9 20])
clear sum_B num den c 
%................mmse
 
 power(1)=6.8;
 for l=2:length(m)
     p=2;
     M=m(l);
     dif=1;
        while dif >.1
             x=10*ones(K,n);
             a=10*ones(K,n);
             p=p-.01;
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
            Rp1=sum(sum(c)/n);
            dif=abs(Rp-Rp1);
        end 
        power(l)=10*log10(p)
 end 
 6
 hold on
plot(m,power,'-')
axis([min(m) max(m) -9 20])
clear sum_B num den c
xlabel('Number of Base Station Antenna')
ylabel('Required Power , Normalized (dB)')
