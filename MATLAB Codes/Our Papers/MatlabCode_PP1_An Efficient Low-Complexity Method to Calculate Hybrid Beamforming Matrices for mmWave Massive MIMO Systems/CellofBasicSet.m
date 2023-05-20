clear
clc
addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' ))
N=11;
[Setb,~,~,~]=GenUniqueSet(N);
Setb=sort(Setb);
S=exp(1j*(2*pi*(0:N-1)/N-pi/(N))).';
Fig=figure;
polarplot(angle(Setb),abs(Setb),'g.')
hold on
polarplot(angle(Setb(1:120)),abs(Setb(1:120)),'b.')
polarplot(angle(Setb(1:80)),abs(Setb(1:80)),'y.')
polarplot(angle(Setb(1:40)),abs(Setb(1:40)),'r.')
% polarplot(angle(S),abs(S),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .9 0])

ax=gca;
ax.RTick = [abs(Setb(40))  abs(Setb(80)) abs(Setb(120)) max(abs(Setb))];
rlim([0 max(abs(Setb))])
ax.ThetaTick = 0:360/N:360;
% ax.ThetaAxis.Color = 'r';
ax.ThetaColor=[.15 .15 .15]; 
ax.RColor=[.15 .15 .15];
ax.LineWidth = 1.5;

ax.Units='centimeters';
ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
Fig.Units='centimeters';
Fig.Position = [10 5 (ax.TightInset(1)+ 2*ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];
