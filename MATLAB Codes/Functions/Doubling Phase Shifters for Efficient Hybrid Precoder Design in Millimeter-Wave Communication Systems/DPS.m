function [WRF,WBB,FRF,FBB]=DPS(H,NRF,Nq)
Ns=NRF;
[U,~,V] = svd(H);
Fopt = V(:,1:Ns);  
Wopt = U(:,1:Ns);
FBB=Fopt(1:NRF,:);
X=Fopt(NRF+1:end,:)/FBB;
FRF=[eye(NRF);X];
Scl=max(max(abs(FRF)))/2;
FRF=FRF/Scl;
FBB=FBB*Scl;
WBB=Wopt(1:Ns,:);
WRF=[eye(NRF);Wopt(NRF+1:end,:)/WBB];
Scl=max(max(abs(WRF)))/2;
WRF=WRF/Scl;
WBB=WBB*Scl;
if (exist('Nq','var'))
    F1=exp(1j*(angle(FRF)+acos(abs(FRF)./2))); 
    F1=Qtz(F1,Nq);
    F2=exp(1j*(angle(FRF)-acos(abs(FRF)./2)));
    F2=Qtz(F2,Nq);
    FRF=F1+F2;
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro'); 
    
    W1=exp(1j*(angle(WRF)+acos(abs(WRF)./2)));
    W1=Qtz(W1,Nq);
    W2=exp(1j*(angle(WRF)-acos(abs(WRF)./2)));
    W2=Qtz(W2,Nq);
    WRF=W1+W2;
%     WBB = sqrt(Ns) * WBB / norm(WRF * WBB,'fro'); 
end




