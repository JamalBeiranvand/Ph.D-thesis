clear
clc
 addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' )) % Add path of the "Functions" folder.
%addpath(('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions\Antenna Selection') )% Add path of the "Functions" folder.

N=6;
[Setb,Sb,Set,~]=GenUniqueSet(N);
S=exp(1j*2*pi*(0:N-1)/N);
max(abs(Set))
Fig=figure;
polarplot(angle(Set(abs(Set)>1.9)),abs(Set(abs(Set)>1.9)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r')
hold on
polarplot(angle(Set(abs(Set)<1.9)),abs(Set(abs(Set)<1.9)),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b')
polarplot(angle(S),abs(S),'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .9 0])
polarplot(0,0,'o','Markersize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','k')
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


Legend1=['Points with one repeatations'];
Legend2=['Points with two repeatations'];
Legend3=['Points with six repeatations'];
Legend4='Point with ten repeatations';
lgd=legend(Legend1,Legend2,Legend3,Legend4);
lgd.Location= 'northwest'; %%
lgd.Box='on';
lgd.LineWidth=.9;
lgd.Interpreter='latex';
lgd.FontSize = 10;
size(Setb)
ax.Units='centimeters';
ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
Fig.Units='centimeters';
Fig.Position = [10 5 (ax.TightInset(1)+ 2*ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))+2];
exportgraphics(gcf,'transparent.eps','ContentType','vector','BackgroundColor','none')

