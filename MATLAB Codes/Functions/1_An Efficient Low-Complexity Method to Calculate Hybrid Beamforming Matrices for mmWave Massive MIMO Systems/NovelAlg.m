function [WRF,WBB,FRF,FBB]=NovelAlg(H,Ns,Fb,Sb)
[U,~,V] = svd(H);
Fopt = V(:,1:Ns);  
Wopt = U(:,1:Ns);
%Ns=size(Fopt,2);
%% Algorithm 
[S,C,~]=BeamSwitch(Fopt,Fb,Sb);
% 7: FRF = SC
    FRF=S*C;
% 8: FBB
    FBB = pinv(FRF) * Fopt;
% 9: Normolize FBB
FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
%% WRF and WBB
[S,C,~]=BeamSwitch(Wopt,Fb,Sb);

% 7: FRF = SC
    WRF=S*C;
% 8: FBB
    WBB = pinv(WRF) * Wopt;

