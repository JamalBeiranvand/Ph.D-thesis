function [WRF,WBB,FRF,FBB]=NovelAlgModified(H,Ns,Fb,Sb)
[U,~,V] = svd(H);
Fopt = V(:,1:Ns);  
Wopt = U(:,1:Ns);
N=size(Sb,2);
b=log2(N);
C=Phases(N,Ns);
Nt=size(Fopt,1);
Ns=size(Fopt,2);
S=zeros(Nt,N*Ns);
%% Algorithm 
% Input
F1=Fopt;
% 1: Normalization
    Scale=max(abs(Fb))/max(max(abs(F1)));
    F1=F1*Scale;
% 2: Compute integer part: M
    OrgAng=angle(F1);
    OrgAng(OrgAng<0)=OrgAng(OrgAng<0)+2*pi;
    M=floor(OrgAng/(2*pi)*N);
% 3: Calculate angle F1
    F=abs(F1).*exp(1j*mod(OrgAng,2*pi/N));
for nt=1:Nt
    L=1;
    for ns=1:Ns
        % 4: Do mapping by using basis set
        Nnn=find(abs(Fb-F(nt,ns))==min(abs(Fb-F(nt,ns))));
        F(nt,ns)=Fb(Nnn(1));
        % 5: Extraction S from code book of the basis set
        Sntns=Sb(Nnn(1),:);
        % 6: Circularly shift S with M
        S(nt,L:ns*N)=circshift(Sntns,M(nt,ns));
        L=ns*N+1;
    end
end
% 7: FRF = SC
    FRF=S*C;
% 8: FBB
    FBB = pinv(FRF) * Fopt;
% 9: Normolize FBB
FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');


%% WRF and WBB
N=size(Sb,2);
% b=log2(N);
Nr=size(Wopt,1);
Ns=size(Wopt,2);
S=zeros(Nr,N*Ns);
%% Algorithm 
% Input
F1=Wopt;
% 1: Normalization
    Scale=max(abs(Fb))/max(max(abs(F1)));
    F1=F1*Scale;
% 2: Compute integer part: M
    OrgAng=angle(F1);
    OrgAng(OrgAng<0)=OrgAng(OrgAng<0)+2*pi;
    M=floor(OrgAng/(2*pi)*N);
% 3: Calculate angle F1
    F=abs(F1).*exp(1j*mod(OrgAng,2*pi/N));
for nt=1:Nr
    L=1;
    for ns=1:Ns
        % 4: Do mapping by using basis set
        Nnn=find(abs(Fb-F(nt,ns))==min(abs(Fb-F(nt,ns))));
        F(nt,ns)=Fb(Nnn(1));
        % 5: Extraction S from code book of the basis set
        Sntns=Sb(Nnn(1),:);
        % 6: Circularly shift S with M
        S(nt,L:ns*N)=circshift(Sntns,M(nt,ns));
        L=ns*N+1;
    end
end
% 7: FRF = SC
    WRF=S*C;
% 8: FBB
    WBB = pinv(WRF) * Wopt;
% 9: Normolize FBB