%......................zf m=100
clear
clc
m=100;
T=196;
K=1:m;
Rp=1:1:58;
for l=1:length(Rp)
    r=Rp(l)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10000;
              dif=1;
              while dif>.000001
                  num=(2^(r*T/(k*(T-t)))-1)*((k+t)*p+1);
                  p1=sqrt(num/(t*(m-k)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          n=find(t_p==min(t_p));
          t_t(j)=min(to(n));
          k_t_p(j)=min(t_p);
          clear n
      end
      n=find(k_t_p==min(k_t_p));
      n=max(n);
      op_t(l)=t_t(n);
      op_k(l)=K(n);
      op_p(l)=min(k_t_p);
      eta(l)=Rp(l)/(.265*op_p(l));
end
plot(Rp,op_t)
  axis([1 60 1 140])
  hold on
  plot(Rp,op_k,'--')
  %..........................mrc m=100
   clear
clc
m=100;
T=196;
K=1:T;
Rp=[1:3:16,18:2:24,26:2:33,35:2:50,52:2:58];
for l=1:length(Rp)
    r=Rp(l)
      for j=1:length(K)
          k=K(j);
          to=k:196;
          t_p=200*ones(1,length(to));
          for i=1:length(to)
              t=to(i);
              p=10000;
              dif=1;
              while dif>10^-6
                  num=(2^(r*T/(k*(T-t)))-1)*((t*(k-1)*p.^2)+(k+t)*p+1);
                  p1=sqrt(num/(t*(m-1)));
                  dif=abs(p-p1);
                  p=p1;
              end 
              t_p(i)=p;
          end 
          n=find(t_p==min(t_p));
          t_t(j)=min(to(n));
          k_t_p(j)=min(t_p);
          clear n
      end
      n=find(k_t_p==min(k_t_p));
      n=max(n);
      op_t(l)=t_t(n);
      op_k(l)=K(n);
      op_p(l)=min(k_t_p);
      eta(l)=Rp(l)/(.265*op_p(l));
      clear t_t k_t_p
end
hold on
plot(Rp,op_t)
axis([1 60 1 140])
plot(Rp,op_k,'r')
xlabel('Spectral-efficincy(bits/s/Hz)')