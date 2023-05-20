function [FRF,FBB]=MUMISO_DownModified(Fopt,Fb,Sb)

[S,C,~]=BeamSwitch(Fopt,Fb,Sb);
% 7: FRF = SC
    FRF=S*C;
% 8: FBB
    FBB = pinv(FRF) * Fopt;
% 9: Normolize FBB
% FBB = sqrt(K) * FBB / norm(FRF * FBB,'fro');
FBB =  FBB / norm(FRF * FBB,'fro');