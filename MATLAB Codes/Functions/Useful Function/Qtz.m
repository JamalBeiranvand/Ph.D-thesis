 % quantize the elements of the X  to the nearest points in the feasible set 
 % input: 
 % X: input matrix (analog precoder matrix or analog combiner matrix)
 % b: the number of bits in the resolution of phase shifters
 
 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en

 function Xq=Qtz(X,Nq) 
 if isnumeric(Nq)
%     Nps=2^b;
    Set=exp(1j*2*pi*(0:Nq-1)/Nq);
    [MatX,MatSet]=meshgrid(X,Set);
    Distance=abs(MatSet-MatX);
    N=(Distance==min(Distance));
    Nn=find((sum(N)~=1));
    for i=1:length(Nn)
        nn=find(N(:,Nn(i))==1);
        N(nn(2:end),Nn(i))=0;
    end
    out=MatSet(N);
    Xq=reshape(out,size(X));
 end
 if ischar(Nq)
     if strcmp(Nq,'inf')
         Xq=X;
     end
 end
 end
 

