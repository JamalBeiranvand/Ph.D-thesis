function [h,power_matrix,A_BS,A_MS]=mmWave_channel(Nr,Nt,L,lamada)
power=sqrt(Nr*Nt/L);
alpha=(randn(1,L)+1i*randn(1,L))/sqrt(2);
power_matrix=power*diag(alpha);
AoA=2*pi*rand(1,L)-pi;
AoD=2*pi*rand(1,L)-pi;
% AoA=(pi/3)*rand(1,L)-pi/6;
% AoD=(pi/3)*rand(1,L)-pi/6;
d=lamada/2;
for l=1:L
    A_BS(:,l)=array_respones(AoD(l),Nt,d,lamada);
    A_MS(:,l)=array_respones(AoA(l),Nr,d,lamada);
end
h=A_MS*power_matrix*A_BS';