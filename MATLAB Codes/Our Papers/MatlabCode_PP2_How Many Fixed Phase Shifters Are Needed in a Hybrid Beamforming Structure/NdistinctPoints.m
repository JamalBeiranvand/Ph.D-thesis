function Nf=NdistinctPoints(N)
P=unique(factor(N));
Xi=N-N/min(P);
if isprime(N)
    A=[];
else
    A=P;
end

Z1=0;
for n=Xi+1:N
    Z1=Z1+nchoosek(N,n);
end

Z2=0;
if  (~isempty(A))
    for k=1:length(A)
        x=A(k);
        Z2=Z2+g(N,Xi,x);
    end
end
C=0;
if N==12 || N==20
    x1=max(A);
%     for z=1:floor(Xi/x1)
%         for M=floor((Xi-z*x1)/x1):Xi-z*x1
%             for n=0:Xi-z*x1-M
%                C=C+(-1)^(z+1)*nchoosek(N/(z*x1),1)*nchoosek(z*x1,M)*nchoosek(N-2*z*x1,n);
%             end
%         end        
%     end
    for M=1:x1
        for n=0:Xi-x1-M
           C=C+nchoosek(N/((floor(M/x1)+1)*x1),1)*nchoosek(x1,M)*nchoosek(N-2*x1,n);
        end
%            C=C-nchoosek(N/((floor(M/x1)+1)*x1),1)*nchoosek(x1,M)*(floor(M/x1)+1)/2*g(N-2*x1,Xi-x1-M,x1);
    end 
     C=C+nchoosek(N/x1,1)*g(N-2*x1,Xi-x1,2);
elseif N==15 || N==21
    x1=max(A);
    x2=min(A);
    for M=1:floor((Xi-x1)/(x2-1))
        for n=0:((Xi-x1)-(x2-1)*M)
            C=C+(-1)^(M+1)*nchoosek(N/x1,1)*nchoosek(x1,M)*nchoosek(N-x1-(x2-1)*M,n);
        end
    end
% elseif N==18
%      x1=max(A);
%     for M=1:x1
%         for n=0:Xi-x1-M
%            C=C+nchoosek(N/(x1),1)*nchoosek(x1,M)*nchoosek(N-2*x1,n);
%         end
% %            C=C-nchoosek(N/((floor(M/x1)+1)*x1),1)*nchoosek(x1,M)*(floor(M/x1)+1)/2*g(N-2*x1,Xi-x1-M,x1);
%     end  
%      C=C+nchoosek(N/x1,1)*g(N-2*x1,Xi-x1,2)
%      
%      V1=g(N-2*x1,5,3);
%      V2=g(N-2*x1,4,3);
%      V3=6*g(N-2*x1,3,3)
%      
%      B=nchoosek(12,0)+nchoosek(12,1)+nchoosek(12,2)+(nchoosek(12,3))
%      V=18*V1+18*V2
%      C=C-V-V3-2
 end
Z=Z1+Z2-C;

N1=0;
N2=0;
Nz=0;
if ~mod(N,2)&& length(A)==2
    X=max(A);
    for r=(X+1)/2:X-1
        N1=N1+nchoosek(N/X,1)*nchoosek(X,r)*(g(N-2*X,3,2)+g(N-2*X,3,3));
        for u=0:Xi-X
           N2=N2+nchoosek(N/X,1)*nchoosek(X,r)*(nchoosek(N-2*X,u));
        end
    end
    if N==12
       Nz=N2-N1-36;
    else
        Nz=N2-N1;
    end
end


Nf=2^N-Z-Nz;
% Nf=2^N-Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function G=g(Nn,XI,Xx)
        G=0;
        if Nn>=XI && XI>=Xx
            for m=1:floor(XI/Xx)
                for t=0:XI-Xx*m
                    G=G+(-1)^(m+1)*nchoosek(Nn/Xx,m)*nchoosek(Nn-Xx*m,t);               
                end
            end
        end
    end
end
