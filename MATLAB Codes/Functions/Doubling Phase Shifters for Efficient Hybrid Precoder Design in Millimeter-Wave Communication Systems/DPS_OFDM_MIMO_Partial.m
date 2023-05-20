function [WBB,WRF,FRF,FBB]=DPS_OFDM_MIMO_Partial(Wopt,Fopt,NRFTx,NRFRx,P)
[Nt,Ns,F]=size(Fopt);
[Nr,~,~]=size(Wopt);
FRF=zeros(Nt,NRFTx);
FOP=[];
    Wk=[];
    for f=1:F
        FOP=[FOP,Fopt(:,:,f)];
        Wk=[Wk,Wopt(:,:,f)];
    end
    WOP(:,:)=Wk;
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
for f=1:F
    SColumn=(f-1)*Ns+1;
    EColumn=f*Ns;
    FBB(:,:,f)=FB(:,SColumn:EColumn);
    FBB(:,:,f)=sqrt(P)*FBB(:,:,f)/norm(FRF*FBB(:,:,f),'fro');
end


%%

    WB=zeros(NRFRx,size(WOP,2));
    for j=1:NRFRx
        Sigyyh=zeros(size(WOP,2));
        for i= (j-1)*(Nr/NRFRx)+1:j*(Nr/NRFRx)
            yi=WOP(i,:).';
            Sigyyh=Sigyyh+yi*yi';
        end
        [V,Lam]=eig(Sigyyh);
        xj=V(:,end);

        for i= (j-1)*(Nr/NRFRx)+1:j*(Nr/NRFRx)
            yi=WOP(i,:).';
            ai=(xj'*yi/(xj'*xj));
            WRF(i,j)=ai;       
        end    
        WB(j,:)= xj.';   
    end
    
    for f=1:F
        SColumn=(f-1)*Ns+1;
        EColumn=f*Ns;
        WBB(:,:,f)=WB(:,SColumn:EColumn);
    end
    
end


