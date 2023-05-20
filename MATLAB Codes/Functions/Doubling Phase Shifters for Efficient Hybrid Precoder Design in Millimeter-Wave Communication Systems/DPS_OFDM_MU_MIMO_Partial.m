function [WBB,WRF,FRF,FBB]=DPS_OFDM_MU_MIMO_Partial(Wopt,Fopt,NRFTx,NRFRx,P)
[Nt,Ns,K,F]=size(Fopt);
[Nr,~,~,~]=size(Wopt);
FRF=zeros(Nt,NRFTx);
FOP=[];
  for k=1:K
        Wk=[];
        for f=1:F
            FOP=[FOP,Fopt(:,:,k,f)];
            Wk=[Wk,Wopt(:,:,k,f)];
        end
        WOP(:,:,k)=Wk;
  end
  FB=zeros(NRFTx,size(FOP,2));
for j=1:NRFTx
    Sigyyh=zeros(size(FOP,2));
    for i= (j-1)*(Nt/NRFTx)+1:j*(Nt/NRFTx)
        yi=FOP(i,:).';
        Sigyyh=Sigyyh+yi*yi';
    end
    [V,Lam]=eig(Sigyyh);
    xj=(V(:,end));
    
    for i= (j-1)*(Nt/NRFTx)+1:j*(Nt/NRFTx)
        yi=FOP(i,:).';
        ai=(xj'*yi/(xj'*xj));
        FRF(i,j)=ai;       
    end    
    FB(j,:)= xj.';   
end
% sum(sum(abs(FOP-FRF*FB).^2))
for k=1:K
    for f=1:F
        SColumn=(k-1)*F*Ns+(f-1)*Ns+1;
        EColumn=(k-1)*F*Ns+f*Ns;
        FBB(:,:,k,f)=FB(:,SColumn:EColumn);
        FBB(:,:,k,f)=sqrt(Ns*P)*FBB(:,:,k,f)/norm(FRF*FBB(:,:,k,f),'fro');
    end
end

%%

for k=1:K
    WB=zeros(NRFRx,size(WOP,2));
    for j=1:NRFRx
        Sigyyh=zeros(size(WOP,2));
        for i= (j-1)*(Nr/NRFRx)+1:j*(Nr/NRFRx)
            yi=WOP(i,:,k).';
            Sigyyh=Sigyyh+yi*yi';
        end
        [V,Lam]=eig(Sigyyh);
        xj=V(:,end);

        for i= (j-1)*(Nr/NRFRx)+1:j*(Nr/NRFRx)
            yi=WOP(i,:,k).';
            ai=(xj'*yi/(xj'*xj));
            WRF(i,j,k)=ai;       
        end    
        WB(j,:)= xj.';   
    end
    
    for f=1:F
        SColumn=(f-1)*Ns+1;
        EColumn=f*Ns;
        WBB(:,:,k,f)=WB(:,SColumn:EColumn);
    end
    
end


