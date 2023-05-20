function [F,Fn,Fc]=SPF(H,SetK)
  [Nt,K]=size(H);
    if Nt<K
        error('size of the channel matrix is not true')
    end
    Nuk=sum(SetK~=0);
    Nu=sum(Nuk);
    Ncom=Nt-Nu;
        %% comput Fn
    Fn=zeros(sum(Nuk),K);
    nt=1;
    for k=1:K    
        for nuk=1:Nuk(k)
            Fn(nt,k)=1;
            nt=nt+1;
        end
    end   
    %% Compute Fc   
    Fc=ones(Ncom,K);
    F=[Fn;Fc];
    %%

    Hc=H(end-Ncom+1:end,:);
    I=eye(K);
   for k=1:K
        Nk=SetK(:,k);
        Nk(Nk==0)=[]; 
        Ik=I(:,k);
        Hnk=H(Nk,:);
        Hk=[Hnk;Hc];       
        fk=conj(Hk)*inv(Hk.'*conj(Hk))*Ik;
        F(F(:,k)==1,k)=fk;       
%         HH(:,:,k)=Hk;
   end
    
   F=F/sqrt(sum(sum(abs(F).^2)));
   Fn=F(1:Nu,:);
   Fc=F(Nu+1:end,:);
end
    
    
    