clear
clc
 addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' )) % Add path of the "Functions" folder.
NPS=9:12;
[Setb1,~,Set1,Sb1]=GenUniqueSet(NPS(1));
[Setb2,~,Set2,Sb2]=GenUniqueSet(NPS(2));
[Setb3,~,Set3,Sb3]=GenUniqueSet(NPS(3));
[Setb4,~,Set4,Sb4]=GenUniqueSet(NPS(4));
AgS1=angle(Set1);AbS1=abs(Set1);
AgS2=angle(Set2);AbS2=abs(Set2);
AgS3=angle(Set3);AbS3=abs(Set3);
AgS4=angle(Set4);AbS4=abs(Set4);
%%
Fig=figure ;
 subplot(2,2,1)
    polarscatter(AgS1,AbS1,'b.')
    hold on
    polarscatter(angle(Setb1),abs(Setb1),'r.')
    ax=gca;
    ax.RTick = [0  4];
%     ax.RTick = [0 max(abs(Setb1)) 4];
%     ax.RTick = [0 max(abs(Setb1))];
%     ax.RTickLabel =[0 max(abs(Setb1))];
     rlim([0 4])
    ax.ThetaTick = 0:360/NPS(1):360;
    ax.ThetaTickLabel ={'0'; num2str(360/NPS(1)); ''; ''; ''; '';'';'';''};
    ax.ThetaAxis.Color = 'r';
    ax.ThetaColor=[.15 .15 .15]; 
    ax.RColor=[.15 .15 .15];
    ax.LineWidth = 1.5;
    ax.Units='centimeters';
    T11=ax.TightInset(1);T12=ax.TightInset(2);
    T13=ax.TightInset(3);T14=ax.TightInset(4);  
    ax.Position = [T11+T13 6 5 5];
    
    Not1=annotation('textbox');
    Not1.String=['$N_{\mathrm{ps}}=$ ',num2str(NPS(1)), ', $\eta_{9}= $',num2str(max(abs(Setb1)))];
    Not1.FontSize = 10;
    Not1.Interpreter='latex';
    Not1.Units='centimeters';
    Not1.Position=[ .8 10.8   3 1];
    Not1.FitBoxToText='on';
    %%
subplot(2,2,2)
    polarscatter(AgS2,AbS2,'b.')
    hold on
    polarscatter(angle(Setb2),abs(Setb2),'r.')
      ax=gca;
      ax.RTick = [0  4];
%     ax.RTick = [0 max(abs(Setb2)) 4];
%     ax.RTick = [0 max(abs(Setb2))];
%     ax.RTickLabel =[0 max(abs(Setb2))];
     rlim([0 4])
    ax.ThetaTick = 0:360/NPS(2):360;
    ax.ThetaTickLabel ={'0'; num2str(360/NPS(2)); ''; ''; ''; '';'';'';'';''};
    ax.ThetaAxis.Color = 'r';
    ax.ThetaColor=[.15 .15 .15]; 
    ax.RColor=[.15 .15 .15];
    ax.LineWidth = 1.5;
    ax.Units='centimeters';
    ax.Position = [6+T11+T13 6 5 5];
    
    Not1=annotation('textbox');
    Not1.String=['$N_{\mathrm{ps}}=$ ',num2str(NPS(2)), ', $\eta_{10}= $',num2str(max(abs(Setb2)))];
    Not1.FontSize = 10;
    Not1.Interpreter='latex';
    Not1.Units='centimeters';
    Not1.Position=[ .8+6 10.8  3 1];
    Not1.FitBoxToText='on';
    %%
subplot(2,2,3)
    polarscatter(AgS3,AbS3,'b.')
    hold on
    polarscatter(angle(Setb3),abs(Setb3),'r.')
    ax=gca;
    ax.RTick = [0  4];
%     ax.RTick = [0 max(abs(Setb3)) 4];
%     ax.RTickLabel =[0 max(abs(Setb3))];
     rlim([0 4])
    ax.ThetaTick = 0:360/NPS(3):360;
    ax.ThetaTickLabel ={'0'; num2str(360/NPS(3)); ''; ''; ''; '';'';'';'';'';''};
% ax.ThetaTick = [0,360/NPS(3)];
%     ax.ThetaTickLabel =[0,360/NPS(3)];
    ax.ThetaAxis.Color = 'r';
    ax.ThetaColor=[.15 .15 .15]; 
    ax.RColor=[.15 .15 .15];
    ax.LineWidth = 1.5;
    ax.Units='centimeters';
    ax.Position = [T11+T13 0 5 5];
    
     Not1=annotation('textbox');
    Not1.String=['$N_{\mathrm{ps}}=$ ',num2str(NPS(3)), ', $\eta_{11}= $',num2str(max(abs(Setb3)))];
    Not1.FontSize = 10;
    Not1.Interpreter='latex';
    Not1.Units='centimeters';
    Not1.Position=[ .8  4.8  3 1];
    Not1.FitBoxToText='on';
    %%
subplot(2,2,4)
    polarscatter(AgS4,AbS4,'b.')
    hold on
    polarscatter(angle(Setb4),abs(Setb4),'r.')
    ax=gca;
    ax.RTick = [0  4];
%     ax.RTick = [0 max(abs(Setb4))];
    ax.RTickLabel =[0  4];
    rlim([0 4])
    ax.ThetaTick = 0:360/NPS(4):360;
    ax.ThetaTickLabel ={'0'; num2str(360/NPS(4)); ''; ''; ''; '';'';'';'';'';'';''};
    ax.ThetaAxis.Color = 'r';
    ax.ThetaColor=[.15 .15 .15]; 
    ax.RColor=[.15 .15 .15];
    ax.LineWidth = 1.5;
    ax.Units='centimeters';
    ax.Position = [6+T11+T13 0 5 5];
    
   Not1=annotation('textbox');
    Not1.String=['$N_{\mathrm{ps}}=$ ',num2str(NPS(4)), ', $\eta_{12}= $',num2str(max(abs(Setb4)))];
    Not1.FontSize = 10;
    Not1.Interpreter='latex';
    Not1.Units='centimeters';
    Not1.Position=[ .8+6 4.8  3 1];
    Not1.FitBoxToText='on';
%%
% Legend1=['Elements of the feasible set for $N=$',num2str(N)];
% Legend2=['Elements of the basic set for $N=$',num2str(N)];
% Legend3='Phase shifters on unit circle';
% lgd=legend(Legend1,Legend2,Legend3);
% lgd.Location= 'northwest'; %%
% lgd.Box='on';
% lgd.LineWidth=.9;
% lgd.Interpreter='latex';
% lgd.FontSize = 10;
% size(Setb)
%  ax.Units='centimeters';
%  ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
Fig.Units='centimeters';
Fig.Position = [10 5 2*(2*T11+2*T13+5)+1.3 12];
exportgraphics(gcf,'transparent.eps','ContentType','vector','BackgroundColor','none')