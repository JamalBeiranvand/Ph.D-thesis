function A=SwitchMatrix(Nt,ComSet,SetK)
    NC=length(ComSet);
    Nuk=sum(SetK~=0);
    K=length(Nuk);
%     Nu=sum(Nuk);
    A=zeros(Nt,Nt);
    i=1;
    for k=1:K    
        for nuk=1:Nuk(k)
            A(SetK(nuk,k),i)=1;
            i=i+1;
        end
    end

    for i=1:NC
         A(ComSet(i),Nt-NC+i)=1;
    end
end