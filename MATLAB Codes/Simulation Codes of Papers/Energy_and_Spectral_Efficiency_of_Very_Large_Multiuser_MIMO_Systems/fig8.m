clear
clc
n=10^3;
P=(.01:.1:500);
Z=randn(n,1);
p=repmat(P,n,1);
z=repmat(Z,1,length(P));
K=1;
to=K:195;
T=196;
m=1;
Rp=1:7;
eta(1)=2.0457;
for l=2:length(Rp)
r=Rp(l);
reference=l
 for i=1:length(to)
         t=to(i);
     num=t*p.^2.*abs(z).^2;
    den=1+p*(1+t);
    x=num./den;
    rp=((T-t)/T)*mean(log2(1+x));
    dif=abs(r-rp);
    m=find(dif==min(dif));
    op_p(i)=min(P(m));
    eta_t(i)=r/(.26*op_p(i));
 end 
 eta(l)=max(eta_t)
end 
semilogy(Rp,eta,':')
 axis([1 90 .1 10^4])
% eta=[ 2.0457  1.0523 0.5626 0.3046 0.1658 0.0904 0.0539]
%..............zf b=.32
 clear
clc
   L=7;
  b=.32;
   l=(L-1)*b+1;
   m=100;
T=196;
K=1:m;
Rp=[1:2:11 11.9780];
for z=1:length(Rp)
    r=Rp(z);
    zf32=Rp(z)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10;
              dif=1;
              while dif>.00001
                  c=(2^(r*T/(k*(T-t)))-1);
                  num=c*(t*k*(l^2+b-(l*b+1))*p^2+l*(k+t)*p+1);
                  p1=sqrt(num/(t*(m-k)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          k_t_p(j)=min(t_p);
          clear to 
      end
      min(k_t_p)
      op_p(z)=min(k_t_p);
      eta(z)=Rp(z)/(.265*op_p(z));
end
hold on
semilogy(Rp,eta,'g')
  axis([1 90 .1 10^4])
  %................................zf b=.11
 clear
clc
   L=7;
  b=.11;
   l=(L-1)*b+1;
   m=100;
T=196;
K=1:m;
Rp=[1:2:11 13:3:32 32.205];
for z=1:length(Rp)
    r=Rp(z);
    zf11=Rp(z)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10;
              dif=1;
              while dif>.00001
                  c=(2^(r*T/(k*(T-t)))-1);
                  num=c*(t*k*(l^2+b-(l*b+1))*p^2+l*(k+t)*p+1);
                  p1=sqrt(num/(t*(m-k)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          k_t_p(j)=min(t_p);
          clear to 
      end
      min(k_t_p)
      op_p(z)=min(k_t_p);
      eta(z)=Rp(z)/(.265*op_p(z));
end
hold on
semilogy(Rp,eta,'r')
  axis([1 90 .1 10^4])
  %.......................zf b=.04
   clear
clc
   L=7;
  b=.04;
   l=(L-1)*b+1;
   m=100;
T=196;
K=1:m;
Rp=[1:3:61 62 62.52];
for z=1:length(Rp)
    r=Rp(z);
    zf04=Rp(z)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10;
              dif=1;
              while dif>.00001
                  c=(2^(r*T/(k*(T-t)))-1);
                  num=c*(t*k*(l^2+b-(l*b+1))*p^2+l*(k+t)*p+1);
                  p1=sqrt(num/(t*(m-k)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          k_t_p(j)=min(t_p);
          clear to 
      end
      min(k_t_p)
      op_p(z)=min(k_t_p);
      eta(z)=Rp(z)/(.265*op_p(z));
end
hold on
semilogy(Rp,eta)
  axis([1 90 .1 10^4])
  %...........................mrc b=.32
   clear
clc
   L=7;
  b=.32;
   l=(L-1)*b+1;
   m=100;
T=196;
K=1:T;
Rp=[1:2:10  10 10.1550];
for z=1:length(Rp)
    r=Rp(z);
    mrc32=Rp(z)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10;
              dif=1;
              while dif>.00001
                  c=(2^(r*T/(k*(T-t)))-1);
                  num=c*((k*l^2-1+b*(l-1)*(m-2))*t*p^2+l*(k+t)*p+1);
                  p1=sqrt(num/(t*(m-1)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          k_t_p(j)=min(t_p);
          clear to 
      end
      min(k_t_p)
      op_p(z)=min(k_t_p);
      eta(z)=Rp(z)/(.265*op_p(z));
end
hold on
semilogy(Rp,eta,'g')
  axis([1 90 .1 10^4])
  %..........................mrc b=.11
   clear
clc
   L=7;
  b=.11;
   l=(L-1)*b+1;
   m=100;
T=196;
K=1:T;
Rp=[1:2:28 28 28.3150];
for z=1:length(Rp)
    r=Rp(z);
     mrc11=Rp(z)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10;
              dif=1;
              while dif>.00001
                  c=(2^(r*T/(k*(T-t)))-1);
                  num=c*((k*l^2-1+b*(l-1)*(m-2))*t*p^2+l*(k+t)*p+1);
                  p1=sqrt(num/(t*(m-1)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          k_t_p(j)=min(t_p);
          clear to 
      end
      min(k_t_p)
      op_p(z)=min(k_t_p);
      eta(z)=Rp(z)/(.265*op_p(z));
end
hold on
semilogy(Rp,eta,'r')
  axis([1 90 .1 10^4])
  %...........................mrc m=.04
   clear
clc
   L=7;
  b=.04;
   l=(L-1)*b+1;
   m=100;
T=196;
K=1:T;
Rp=[1:2:44 44 44.506];
for z=1:length(Rp)
    r=Rp(z);
     mrc04=Rp(z)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10;
              dif=1;
              while dif>.00001
                  c=(2^(r*T/(k*(T-t)))-1);
                  num=c*((k*l^2-1+b*(l-1)*(m-2))*t*p^2+l*(k+t)*p+1);
                  p1=sqrt(num/(t*(m-1)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          k_t_p(j)=min(t_p);
          clear to 
      end
      min(k_t_p)
      op_p(z)=min(k_t_p);
      eta(z)=Rp(z)/(.265*op_p(z));
end
hold on
semilogy(Rp,eta)
  axis([1 90 .1 10^4])
  %.............................
   clear
clc
P=-20:10:20;
rp=1:90;
for z=1:length(P)
  p=10^(P(z)/10);
   eta=rp/(.265*p);
   hold on
   semilogy(rp,eta,':')
end
  axis([1 90 .1 10^4])
 xlabel('spectral-efficincy(bits/s/Hz)')
 ylabel('Relative Energy-efficincy (bits/J)')