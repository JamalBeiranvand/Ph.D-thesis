clear
clc
 addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' )) % Add path of the "Functions" folder.
%addpath(('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions\Antenna Selection') )% Add path of the "Functions" folder.

N=10;
[Setb,Sb,Set,~]=GenUniqueSet(N);

Set=exp(1j*2*pi*(0:N-1)/N);
S=exp(1j*2*pi*(0:3)/N);
max(abs(Set))
Fig=figure;
subplot(1,2,1)
%polarplot(angle(Set(abs(Set)>1.1)),abs(Set(abs(Set)>1.1)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r')
%hold on
polarplot(angle(Set(abs(Set)<1.1)),abs(Set(abs(Set)<1.1)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0])
hold on
polarplot(angle(S),abs(S),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .9 0])
polarplot(angle(sum(S)),abs(sum(S)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r')
% Phi=0:.1:2*pi;
% polarplot(Phi,ones(length(Phi),1),'Color',[.15 .15 .15] )
ax=gca;
ax.RTick = [0  1 1.732 max(abs(Setb))];
 ax.RTickLabel =[];
rlim([0 2])
ax.ThetaTick = 0:360/N:360;
ax.ThetaTickLabel={'0°'; '60°'; '120°'; '180°'; '240°'; '300°'; '360°'};
ax.ThetaAxis.Color = 'r';
ax.ThetaColor=[.15 .15 .15]; 
ax.RColor=[.15 .15 .15];
ax.LineWidth = 1;
% ax.TightInset='top';


Legend1=['Inactive PSs'];
Legend2=['Active PSs'];
Legend3=['Generated coefficient'];
Legend4='The zero point';
lgd=legend(Legend1,Legend2,Legend3,Legend4);
lgd.Location= 'northwest'; %%
lgd.Box='on';
lgd.LineWidth=.9;
lgd.Interpreter='latex';
lgd.FontSize = 10;
size(Setb)
% ax.Units='centimeters';
% ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
% Fig.Units='centimeters';
% Fig.Position = [10 5 (ax.TightInset(1)+ 2*ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))+2];
%
N=6;
Set=exp(1j*2*pi*(0:N-1)/N);
S=exp(1j*2*pi*(1:2)/N);
subplot(1,2,2)
%polarplot(angle(Set(abs(Set)>1.1)),abs(Set(abs(Set)>1.1)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r')
%hold on
polarplot(angle(Set(abs(Set)<1.1)),abs(Set(abs(Set)<1.1)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0])
hold on
polarplot(angle(S),abs(S),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .9 0])
polarplot(angle(sum(S)),abs(sum(S)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r')
% Phi=0:.1:2*pi;
% polarplot(Phi,ones(length(Phi),1),'Color',[.15 .15 .15] )
ax=gca;
ax.RTick = [0  1 1.732 max(abs(Setb))];
 ax.RTickLabel =[];
rlim([0 2])
ax.ThetaTick = 0:360/N:360;
ax.ThetaTickLabel={'0°'; '60°'; '120°'; '180°'; '240°'; '300°'; '360°'};
ax.ThetaAxis.Color = 'r';
ax.ThetaColor=[.15 .15 .15]; 
ax.RColor=[.15 .15 .15];
ax.LineWidth = 1;
% ax.TightInset='top';


Legend1=['Inactive PSs'];
Legend2=['Active PSs'];
Legend3=['Generated coefficient'];
Legend4='The zero point';
lgd=legend(Legend1,Legend2,Legend3,Legend4);
lgd.Location= 'northwest'; %%
lgd.Box='on';
lgd.LineWidth=.9;
lgd.Interpreter='latex';
lgd.FontSize = 10;
size(Setb)