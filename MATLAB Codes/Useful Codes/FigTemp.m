% this is a good template to plot figures! 
% Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
% google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en

% prepare data to plot
Xaxis=0:10;
Power=[2,1.7,1.5];
Data1=Xaxis.^2;
Data2= Xaxis.^1.7;
Data3=Xaxis.^1.5;
%% plot 
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
Legend1=['Data1 $y=x^{',num2str(Power(1)),'}$'];
Legend2=['Data1 $y=x^{',num2str(Power(2)),'}$'];
Legend3=['Data1 $y=x^{',num2str(Power(3)),'}$'];

set(P1,'LineStyle',':');           % '-' , '--' , ':' , '-.', 'non'
set(P1,'Marker','hexagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P1,'LineWidth',2);
set(P1,'Color','b');
set(P1,'MarkerSize',7);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

set(P2,'LineStyle','--');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color','r');
set(P2,'MarkerSize',7);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

set(P3,'LineStyle','-.');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','pentagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color','k');
set(P3,'MarkerSize',7);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

lgd=legend(Legend1,Legend2,Legend3);

% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
% ax.YLim = [5 40];
ax.XTick = Xaxis;
ax.XTickLabel=Xaxis;
ax.XLabel.String='$\alpha$ label';
ax.YLabel.String='y label';
ax.XGrid='on';
ax.YGrid='on';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'on';
ax.TickLabelInterpreter='latex';
ax.YLabel.Interpreter='latex';
ax.XLabel.Interpreter='latex';
ax.FontSize = 11;
ax.Units='centimeters';
ax.LineWidth = .9;
ax.Box='on';
% Legend Properties %
lgd.Location= 'northwest'; %%
lgd.Box='on';
lgd.LineWidth=.9;
lgd.Interpreter='latex';
lgd.FontSize = 10;

% Figure Properties
Fig.Units='centimeters';
ax.Position = [ax.TightInset(1) ax.TightInset(2) 11 8];
Fig.Position = [5 2 (ax.TightInset(1)+ ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];

