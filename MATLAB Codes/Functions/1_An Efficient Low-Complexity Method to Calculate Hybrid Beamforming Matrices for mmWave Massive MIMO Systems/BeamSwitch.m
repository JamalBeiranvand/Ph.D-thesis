function [S,C,Scale]=BeamSwitch(Fopt,Fb,Sb)

[Nt,NRF]=size(Fopt);
N=size(Sb,2);
C=Phases(N,NRF);
S=zeros(Nt,N*NRF);
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
    for k=1:NRF
        % 4: Do mapping by using basis set
        Nnn=find(abs(Fb-F(nt,k))==min(abs(Fb-F(nt,k))));
        F(nt,k)=Fb(Nnn(1));
        % 5: Extraction S from code book of the basis set
        Sntns=Sb(Nnn(1),:);
        % 6: Circularly shift S with M
        S(nt,L:k*N)=circshift(Sntns,M(nt,k));
        L=k*N+1;
    end
end
% 7: FRF = SC
%     FRF=S*C;