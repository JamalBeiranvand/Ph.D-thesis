% clear
% clc
% n=10^3;
% P=(.01:.1:500);
% Z=randn(n,1);
% p=repmat(P,n,1);
% z=repmat(Z,1,length(P));
% K=1;
% to=K:195;
% T=196;
% m=1;
% Rp=1:7;
% eta(1)=2.0457;
% for l=2:length(Rp)
% r=Rp(l);
% reference=l
%  for i=1:length(to)
%          t=to(i);
%      num=t*p.^2.*abs(z).^2;
%     den=1+p*(1+t);
%     x=num./den;
%     rp=((T-t)/T)*mean(log2(1+x));
%     dif=abs(r-rp);
%     m=find(dif==min(dif));
%     op_p(i)=min(P(m));
%     eta_t(i)=r/(.26*op_p(i));
%  end 
%  eta(l)=max(eta_t)
% end 
% hold on
% semilogy(Rp,eta,':')
%  axis([1 90 .1 10^4])
% 
% %......................................................
% clear
% clc
% p=.01:.01:1000;
% K=1;
% to=K:195;
% T=196;
% m=100;
% Rp=1:15;
% for l=1:length(Rp)
%     l
%  for i=1:length(to)
%          t=to(i);
%          num=t*p*(m-1);
%          den=(t*p+1)*(K-1)+(t+1)+1./p;
%          x=num./den;
%          c=log2(1+x);
%          Rp1=((T-t)*K/T)*c;
%          dif=abs(Rp(l)-Rp1);
%          n=find(dif==min(dif));
%          power(i)=p(n);
%  end 
%  eta(l)=Rp(l)/(.26*min(power));
% end 
% hold on
% semilogy(Rp,eta,'--')
%  axis([1 90 .1 10^4])
% % %...................zf m=100
%  clear
% clc
% m=100;
% T=196;
% K=1:m;
% Rp=1:1:90;
% for l=1:length(Rp)
%     r=Rp(l);
%     z_100=l
%       for j=1:length(K)
%           k=K(j);
%           to=k:196;
%           t_p=200*ones(1,length(to));
%           for i=1:length(to)
%               t=to(i);
%               p=10000;
%               dif=1;
%               while dif>.000001
%                   num=(2^(r*T/(k*(T-t)))-1)*((k+t)*p+1);
%                   p1=sqrt(num/(t*(m-k)));
%                   dif=abs(p-p1);
%                   p=p1;
%               end 
%               t_p(i)=p;
%           end 
%           k_t_p(j)=min(t_p);
%       end
%       min(k_t_p)
%       op_p(l)=min(k_t_p);
%       eta(l)=Rp(l)/(.265*op_p(l));
% end
% hold on
% semilogy(Rp,eta,'r')
%   axis([1 90 .1 10^4])
%   %................zf m=50
%    clear
% clc
% m=50;
% T=196;
% K=1:m;
% Rp=1:1:90;
% for l=1:length(Rp)
%     r=Rp(l);
%     z_50=l
%       for j=1:length(K)
%           k=K(j);
%           to=k:196;
%           t_p=200*ones(1,length(to));
%           for i=1:length(to)
%               t=to(i);
%               p=10000;
%               dif=1;
%               while dif>.0001
%                   num=(2^(r*T/(k*(T-t)))-1)*((k+t)*p+1);
%                   p1=sqrt(num/(t*(m-k)));
%                   dif=abs(p-p1);
%                   p=p1;
%               end 
%               t_p(i)=p;
%           end 
%           k_t_p(j)=min(t_p);
%       end
%       min(k_t_p)
%       op_p(l)=min(k_t_p);
%       eta(l)=Rp(l)/(.265*op_p(l));
% end
% hold on
% semilogy(Rp,eta,'r')
%   axis([1 90 .1 10^4])
%   %....................................mrc m=100. ........
%  clear
% clc
% m=100;
% T=196;
% K=1:T;
% Rp=[1:16,18:24,26:33,35:50,52:60];
% for l=1:length(Rp)
%     r=Rp(l)
%       for j=1:length(K)
%           k=K(j);
%           to=k:196;
%           t_p=200*ones(1,length(to));
%           for i=1:length(to)
%               t=to(i);
%               p=10000;
%               dif=1;
%               while dif>10^-4
%                   num=(2^(r*T/(k*(T-t)))-1)*((t*(k-1)*p.^2)+(k+t)*p+1);
%                   p1=sqrt(num/(t*(m-1)));
%                   dif=abs(p-p1);
%                   p=p1;
%               end 
%               t_p(i)=p;
%           end 
%           n=find(t_p==min(t_p));
%           k_t_p(j)=min(t_p);
%       end
%       op_p(l)=min(k_t_p);
%       eta(l)=Rp(l)/(.265*op_p(l));
%       clear t_t k_t_p
% end
% hold on
% semilogy(Rp,eta,'g')
% axis([1 90 .1 10^4])
% %.......................mrc m=50........
% clear
% clc
% m=50;
% T=196;
% K=1:T;
% Rp=[1:5 7:11 13:17 19:24 26 28:29 31:35 37:40];
% for l=1:length(Rp)
%     r=Rp(l)
%       for j=1:length(K)
%           k=K(j);
%           to=k:196;
%           t_p=200*ones(1,length(to));
%           for i=1:length(to)
%               t=to(i);
%               p=10000;
%               dif=1;
%               while dif>10^-4
%                   num=(2^(r*T/(k*(T-t)))-1)*((t*(k-1)*p.^2)+(k+t)*p+1);
%                   p1=sqrt(num/(t*(m-1)));
%                   dif=abs(p-p1);
%                   p=p1;
%               end 
%               t_p(i)=p;
%           end 
%           n=find(t_p==min(t_p));
%           k_t_p(j)=min(t_p);
%       end
%       op_p(l)=min(k_t_p);
%       eta(l)=Rp(l)/(.265*op_p(l));
%       clear t_t k_t_p
% end
% hold on
% semilogy(Rp,eta,'g')
% axis([1 90 .1 10^4])
% .................................
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