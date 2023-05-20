clear
clc
N=6
P=unique(factor(N));
Xi=N-N/min(P)+1;
if isprime(N)
    A=[];
else
    A=P;
end
Z=zeros(N+1,1);
C=zeros(N+1,1);
for n=0:N
    
    if (isempty(A)) && (0<=n) && (n<Xi)
        Z(n+1)=0;
    end
    
    if (~isempty(A)) && (0<=n) && (n<min(A))
        Z(n+1)=0;
    end
    
    if (~isempty(A)) && (min(A)<=n) && (n<Xi)
        for i=1:length(A)
            for m=1: floor(n/A(i))
                Z(n+1)=Z(n+1)+(-1)^(m+1)*nchoosek(N/A(i),m)*nchoosek(N-m*A(i),n-m*A(i));
            end
        end
    end
    
    if (Xi<=n)
        Z(n+1)=nchoosek(N,n);
    end   
    %% C
    if (N==15) && (sum(A)-1<=n) && (n<Xi)
       for m=1: floor((n-5)/2)
           C(n+1)=C(n+1)+(-1)^(m+1)*nchoosek(10/2,m)*nchoosek(10-m*2,(n-5)-m*2);
       end
        C(n+1)=nchoosek(3,1)*C(n+1);       
    end
    
     if (N==12) && (sum(A)-1<=n) && (n<Xi)
       for m=1: (n-3)
           C(n+1)=C(n+1)+nchoosek(4/(floor(m/3)+1),1)*nchoosek(3,m)*nchoosek(6,n-3-m);
       end
       for m=1: floor((n-3)/2)
             C(n+1)=C(n+1)+nchoosek(4,1)*(-1)^(m+1)*nchoosek(6/2,m)*nchoosek(6-m*2,(n-3)-m*2);
       end      
    end
    Z(n+1)=Z(n+1)-C(n+1);
end
NZ=zeros(N+1,1);
for n=0:N
    if ((N==6) || (N==10) ||(N==14)) && (((max(A)+1)/2)<=n) && (n<Xi-1)
        NZ(n+1)=nchoosek(N/max(A),1)*nchoosek(max(A),n);
    end
    
    if (N==12)
        
        if n==2
            NZ(n+1)=nchoosek(4,1)*nchoosek(3,2)*nchoosek(6,n-2);  
        end
         
         if n==3
            NZ(n+1)=nchoosek(4,1)*nchoosek(3,2)*nchoosek(6,n-2);  
         end
        
        if n==4
            NZ(n+1)=nchoosek(4,1)*nchoosek(3,2)*nchoosek(6,n-2);
           for m=1: floor((n-2)/2)
               NZ(n+1)=NZ(n+1)-12*(-1)^(m+1)*nchoosek(6/2,m)*nchoosek(6-m*2,(n-2)-m*2);
           end  
           NZ(n+1)=NZ(n+1)-nchoosek(2,1)*nchoosek(3,2)*nchoosek(2,1)*nchoosek(3,2);
        end
        
        if n==5
             NZ(n+1)=nchoosek(4,1)*nchoosek(3,2)*nchoosek(6,n-2);
           for m=1: floor((n-2)/2)
               NZ(n+1)=NZ(n+1)-12*(-1)^(m+1)*nchoosek(6/2,m)*nchoosek(6-m*2,(n-2)-m*2);
           end 
           
          for m=1: floor((n-2)/3)
               NZ(n+1)=NZ(n+1)-12*(-1)^(m+1)*nchoosek(6/3,m)*nchoosek(6-m*3,(n-2)-m*3);
          end 
        end
    end
end
S=zeros(N+1,1);
for n=0:N
    S(n+1)=nchoosek(N,n);
end
U=(S-Z-NZ)
Upoints=sum(U)