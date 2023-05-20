% All functions, that are used in this code, are described on the "MATLAb Functions" document.
 % you can see more detail in the functions codes. 
 
 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 14/05/2020
%%
clear;
clc
warning off
%% User pannel Parameters (you may have to change some values in the "Plot" section)
addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions')) % Add path of the "Functions" folder.
K=5;                           % number of users
Ncom=15;
Nps=11;
% b = 3;                    % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 2:2:16;                       % Number of the transmitter antennas on z axis
Txy =16;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 1;                       % Number of the receiver antennas on y axis
% Ns  = 1;                       % Number of data streams
% NrfTx = 9;                     % Number of RF chains at BS
Ncl =5;                       % Number of Channel clusters(Scatters)
Nray= 2;                       % Number of rays in each cluster
AngSpread = 10;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
TxSector=[60 360];
RxSector=[20 360];
SNR_dB = 10;             % Signal to noise ratio in dB
realization = 100;              % Iteration of simulation
%% Initialize Parameters
%Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
%Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);   
%% Allocate memory space for matrices
NT=zeros(size(Txz));
SumR_Fixed=zeros(length(NT),realization);
R_ZF=zeros(length(NT),realization);
R_MRC=zeros(length(NT),realization);
SumR_Dynamic=zeros(length(NT),realization);
DynamicMargin=zeros(length(NT),realization);

%%
[Setb,Sb]=GenUniqueSet(Nps);
for r=1:length(Txz)
    Tx  = [Txz(r),Txy]; 
    NT(r)=Txz(r)*Txy;
    Nt=NT(r);
    Nu=Nt-Ncom;
    disp(r)
    Snr=SNR;
    for reali=1:realization
        % generate channel
        [H,At,Ar,~]=MU_Channel(Tx,K,Ncl,'Nray',Nray,'AngSpread',AngSpread);
        %H=sqrt(1/2)*(randn(K,Nt)+1j*randn(K,Nt));
        G= diag(ones(K,1))*sqrt(Snr);  
        % Full Digital ZF beamforming
        Fzf=ZF_Down(H);
        [FRF_ZF,~]=MUMISO_DownModified(Fzf,Setb,Sb);
        FRF_ZF=FRF_ZF/sqrt(sum(sum(abs(FRF_ZF).^2)));        
        R_ZF(r,reali)=SumRate(H,FRF_ZF,G);
        
%         % Full Digital MRC beamforming
%         Fmrc=MRC_Down(H);
%         R_MRC(r,reali)=SumRate(H,Fmrc,G);        
        %
        H=H.';       
          %% Fixed
        [KSets_Fixed,~]=NonSharedFixed(Nu,K);
        F_Fixed=SPF(H,KSets_Fixed);        
        SumR_Fixed(r,reali)=SumRate(H.',F_Fixed,G);
        [S,C]=BeamSwitch(F_Fixed,Setb,Sb);
        FRF_FSP=S*C;
        FRF_FSP=FRF_FSP/sqrt(sum(sum(abs(FRF_FSP).^2)));
        SumR_Fixed(r,reali)=SumRate(H.',FRF_FSP,G);       
    end
end
Data1= sum(R_ZF,2)/realization;
%Data2= sum(DynamicMargin,2)/realization;
% Data3= sum(SumR_Dynamic,2)/realization;
Data4= sum((SumR_Fixed),2)./realization;
%% plot 

Xaxis= NT;
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
% P2=plot(Xaxis,Data2);
% P3=plot(Xaxis,Data3);
P4=plot(Xaxis,Data4);
% P5=plot(Xaxis,Data5);


Legend1='Fully-connected Structure';
% Legend2=['Dynamic with $Margin=$',num2str(Margin)];
% Legend3='Dynamic Structure';
Legend4='Fixed Structure';

% P1 properties
set(P1,'LineStyle','-');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','*');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color','b');
set(P1,'MarkerSize',7);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

% P2 properties
% set(P2,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P2,'Marker','*');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P2,'LineWidth',2);
% set(P2,'Color',[0 0.5 0.6]);
% set(P2,'MarkerSize',7);
% P2_group= hggroup;set(P2,'Parent',P2_group);
% set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% % P2 properties
% set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P3,'Marker','*');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P3,'LineWidth',2);
% set(P3,'Color',[0 .6 0]);
% set(P3,'MarkerSize',7);
% P3_group= hggroup;set(P3,'Parent',P3_group);
% set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% P2 properties
set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P4,'Marker','s');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P4,'LineWidth',2);
set(P4,'Color',[.8 0 0]);
set(P4,'MarkerSize',11);
P4_group= hggroup;set(P4,'Parent',P4_group);
set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% % P2 properties
% set(P5,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P5,'Marker','x');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P5,'LineWidth',2);
% set(P5,'Color',[0 .5 0]);
% set(P5,'MarkerSize',7);
% P5_group= hggroup;set(P5,'Parent',P5_group);
% set(get(get(P5_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');


lgd=legend(Legend1,Legend4);
% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
 ax.YLim = [10 max(Data1)];
ax.XTick = Xaxis;
% ax.YTick = 0:5:50;
ax.XTickLabel=Xaxis;
ax.XLabel.String='$N_t$';
ax.YLabel.String='Sum Rate (bits/sec/Hz)';
ax.XGrid='on';
ax.YGrid='on';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'on';
ax.YLabel.Interpreter='latex';
ax.XLabel.Interpreter='latex';
ax.FontSize = 11;
ax.LineWidth = .9;
ax.Box='on';

% Legend Properties %
lgd.Location= 'northwest'; %%
lgd.Box='on';
lgd.LineWidth=.9;
lgd.Interpreter='latex';
lgd.FontSize = 10;
% 
% Figure Properties
ax.Units='centimeters';
ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
Fig.Units='centimeters';
Fig.Position = [10 5 (ax.TightInset(1)+ 2*ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];

%%  Creat notation 
Not1=annotation('textbox');
X_points=[ .7*(diff(ax.XLim)) .9*(diff(ax.XLim))]+ax.XLim(1);
Y_points=[ .1*(diff(ax.YLim)) .3*(diff(ax.YLim))]+ax.YLim(1);
str1=['$SNR=$ ',num2str(SNR_dB)];
str2=['$K=$',num2str(K)];
str3=['$N_{com}=$ ',num2str(Ncom)];
Not1.String={str1,str2,str3};
Not1.FontSize = 10;
Not1.Interpreter='latex';
Not1.Units='centimeters';
Dx=ax.Position(3)/(diff(ax.XLim));
Dy=ax.Position(4)/(diff(ax.YLim));
x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
Not1.Position=[x_begin y_begin  x_End y_End];
% % 
% Ell=annotation('ellipse');
% X_points=[19 21];
% Y_points=[25 55];
% Ell.Units='centimeters';
% Dx=ax.Position(3)/(diff(ax.XLim));
% Dy=ax.Position(4)/(diff(ax.YLim));
% x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% Ell.Position=[x_begin y_begin  x_End y_End];
% 
% Textarrow1=annotation('textarrow');
% X_points=[8 8];
% Y_points=[30 62];
% str1=['$b_{ps}^{t}=b_{ps}^{r}=\infty$'];
% Textarrow1.String={str1};
% Textarrow1.Interpreter='latex';
% Textarrow1.FontSize = 10;
% Textarrow1.Units='centimeters';
% Dx=ax.Position(3)/(diff(ax.XLim));
% Dy=ax.Position(4)/(diff(ax.YLim));
% x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% Textarrow1.Position=[x_begin y_begin  x_End y_End];
% %%
% Ell2=annotation('ellipse');
% X_points=[19 21];
% Y_points=[90 125];
% Ell2.Units='centimeters';
% Dx=ax.Position(3)/(diff(ax.XLim));
% Dy=ax.Position(4)/(diff(ax.YLim));
% x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% Ell2.Position=[x_begin y_begin  x_End y_End];
% 
% Textarrow2=annotation('textarrow');
% X_points=[16 16];
% Y_points=[30 81];
% str1=['$b_{ps}^{t}=b_{ps}^{r}=$',num2str(B_Vec{2})];
% Textarrow2.String={str1};
% Textarrow2.Interpreter='latex';
% Textarrow2.FontSize = 10;
% Textarrow2.Units='centimeters';
% Dx=ax.Position(3)/(diff(ax.XLim));
% Dy=ax.Position(4)/(diff(ax.YLim));
% x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% Textarrow2.Position=[x_begin y_begin  x_End y_End];
% %%
% Ell3=annotation('ellipse');
% X_points=[19 21];
% Y_points=[140 185];
% Ell3.Units='centimeters';
% Dx=ax.Position(3)/(diff(ax.XLim));
% Dy=ax.Position(4)/(diff(ax.YLim));
% x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% Ell3.Position=[x_begin y_begin  x_End y_End];
% 
% Textarrow3=annotation('textarrow');
% X_points=[24 24];
% Y_points=[30 90];
% str1=['$b_{ps}^{t}=b_{ps}^{r}=$',num2str(B_Vec{3})];
% Textarrow3.String={str1};
% Textarrow3.Interpreter='latex';
% Textarrow3.FontSize = 10;
% Textarrow3.Units='centimeters';
% Dx=ax.Position(3)/(diff(ax.XLim));
% Dy=ax.Position(4)/(diff(ax.YLim));
% x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% Textarrow3.Position=[x_begin y_begin  x_End y_End];