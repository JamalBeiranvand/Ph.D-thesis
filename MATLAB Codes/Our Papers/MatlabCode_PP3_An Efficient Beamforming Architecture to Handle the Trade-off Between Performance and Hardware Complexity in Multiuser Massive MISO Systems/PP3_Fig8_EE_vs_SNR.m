% Figure 8
% 
% This m-file is written according  to  the paper entitled :
%  "An Efficient Beamforming Architecture to Handle the Trade-off Between Performance and Hardware Complexity in Multiuser Massive MISO Systems" 
% you can see the paper on this link: https://ieeexplore.ieee.org/document/9991158

% Authors: 
% 1- Jamal Beiranvand 
%     Website: 
%     e-mail: Jamalbeiranvand@gmail.com
%     google scholar:  https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
% 2- Vahid Meghdadi
%     website: https://www.unilim.fr/pages_perso/vahid/ 
%     e-mail: meghdadi@ensil.unilim.fr
%     google scholar:  https://scholar.google.com/citations?user=_HYFga8AAAAJ&hl=en
% 
% Required Functions to run
% All functions used in this code are in the "Functions" folder, described in the "MATLAB Functions" document.
% Note: change  the path of Functions folder according its path in your PC.
% 
% Date : 16/12/2022
% By Jamal Beiranvand
%%
clear;
clc
warning off
% User pannel Parameters (you may have to change some values in the "Plot" section)
 addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' ))% Add path of the "Functions" folder.

K=10; % number of users
Ng=2;
Nps=11;
% b = 3;                    % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 16;                       % Number of the transmitter antennas on z axis
Txy =16;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 1;                       % Number of the receiver antennas on y axis
% Ns  = 1;                       % Number of data streams
% NrfTx = 9;                     % Number of RF chains at BS
Ncl =10;                       % Number of Channel clusters(Scatters)
Nray= 5;                       % Number of rays in each cluster
AngSpread = 10;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
TxSector=[120 60];
RxSector=[20 360];
SNR_dB = 0:2:20;             % Signal to noise ratio in dB
realization = 60;              % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);  
Ncom=floor((Nt*K/Ng-Nt)/(K-1));
%% Allocate memory space for matrices
SumR_Fixed=zeros(length(SNR_dB),realization);
SumR_ZF=zeros(length(SNR_dB),realization);
SumR_MRC=zeros(length(SNR_dB),realization);
R_DPS=zeros(length(SNR_dB),realization);
R_ShrQ=zeros(length(SNR_dB),realization);
SumR_Dynamic=zeros(length(SNR_dB),realization);
SumR_DynamicMargin=zeros(length(SNR_dB),realization);
%%
[Setb,Sb]=GenUniqueSet(Nps);
Nq=2^3;
for r=1:length(SNR_dB)
    disp(r)
    Nu=Nt-Ncom;
    NSW=Nu+Ncom*(K);
    EFix=SNR(r)+(K*300+NSW*Nps*5+Nps*K*5)*10e-3;
    ED=SNR(r)+(K*300+NSW*Nps*5+Nt*5+Nps*K*5)*10e-3;
    Eful=SNR(r)+(200+K*300+Nt*K*Nps*5+Nps*K*5)*10e-3;
    EG=SNR(r)+(200+K*300+Nt/Ng*K*Nps*5+Nps*K*5)*10e-3;
    EDPS=SNR(r)+(200+K*300+2*Nt*K*50)*10e-3;
    EQPS=SNR(r)+(200+(K+1)*300+Nt*(K+1)*30)*10e-3;
    Snr=SNR(r);
    for reali=1:realization
        % generate channel
        [H,At,Ar,~]=MU_Channel(Tx,K,Ncl,'Nray',Nray,'AngSpread',AngSpread,'TxSector',TxSector);
%         H=sqrt(1/2)*(randn(K,Nt)+1j*randn(K,Nt));
         G= diag(ones(K,1))*sqrt(Snr);   
        % Full Digital ZF beamforming
        Fzf=ZF_Down(H);        
        SumR_ZF(r,reali)=SumRate(H,Fzf,G);
        
        R_DPS(r,reali)=SumRate(H,Fzf,G);
            
        [VRF1,VD1]= Sohrabi_Alg3(H,K+1,Snr);
        R_ShrQ(r,reali)=SumRate(H,VRF1*VD1,eye(K));
        
%         %%
         H=H.';
           %% Fixed
        [KSets_Fixed,~]=NonSharedFixed(Nu,K);
        F_Fixed=SPF(H,KSets_Fixed);       
        SumR_Fixed(r,reali)=SumRate(H.',F_Fixed,G);
        %% Dynamic
        [ComSet,~]=ComAntSelection_Greedy(H,Ncom); 
        KSets=nComAntSelection_Greedy(H,ComSet);
        A=SwitchMatrix(Nt,ComSet,KSets);
        [F_D,HD,~,~]=SPD(H,A,KSets);
        SumR_Dynamic(r,reali)=SumRate(HD.',F_D,G);
        %% FPS Group 
        [FRF_FPSG,FBB_FPSG]=FPS_AltMin_Group_MU(H.',K,K,Nps,Ng);
        SumR_FPG(r,reali)=SumRate(H.',FRF_FPSG*FBB_FPSG,G);
    end
end

%% plot
Data1= sum(SumR_ZF,2)/realization;
Data2=sum(R_DPS,2)/realization;
Data3=sum(R_ShrQ,2)/realization;
Data4= sum(SumR_Dynamic,2)/realization;
Data5= sum((SumR_Fixed),2)/realization;
Data6= sum(SumR_FPG,2)/realization;
% plot 
%%
Xaxis= SNR_dB;
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
P4=plot(Xaxis,Data4);
P5=plot(Xaxis,Data5);
P6=plot(Xaxis,Data6);




Legend6='FC Architecture, DPS [28]';
Legend5= 'FC Architecture [21]';
Legend3='FC Architecture, for 3-bit PSs [16]';
Legend1='Dynamic PFC Architecture Implemented by the Proposed DNN';
Legend2='PFC Architecture';
Legend4=['GC Architecture, for $g=$',num2str(Ng),' [18]'];

% Legend3='Fully Digital ZF';

% Legend4='Proposed Hybrid BF Alg for $b=1$, $N_{RF}=N_s+3 $';
% Legend5=' Quantized-Hybrid BF in [31] (OMP)';




% 

% P2 properties
set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P4,'Marker','pentagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P4,'LineWidth',2);
set(P4,'Color',[.0 0.9 0.0]);
set(P4,'MarkerSize',7);
P4_group= hggroup;set(P4,'Parent',P4_group);
set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P2 properties
set(P5,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P5,'Marker','o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P5,'LineWidth',2);
set(P5,'Color',[0.0 .5 0.0]);
set(P5,'MarkerSize',7);
P5_group= hggroup;set(P5,'Parent',P5_group);
set(get(get(P5_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P2 properties
set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','*');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color',[0.5 .0 0]);
set(P3,'MarkerSize',5);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 

% P2 properties
set(P6,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P6,'Marker','square');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P6,'LineWidth',2);
set(P6,'Color',[0.8 .0 0.0]);
set(P6,'MarkerSize',7);
P6_group= hggroup;set(P6,'Parent',P6_group);
set(get(get(P6_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P2 properties
set(P2,'LineStyle','--');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','hexagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color','b');
set(P2,'MarkerSize',7);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P1 properties
set(P1,'LineStyle','-');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','none');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color','r');
set(P1,'MarkerSize',7);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];




lgd=legend(Legend1,Legend2,Legend3,Legend4,Legend5,Legend6);
% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) 10];
ax.YLim = [0 74];
ax.XTick = Xaxis;
% ax.YTick = 0:5:50;
ax.XTickLabel=Xaxis;
ax.XLabel.String='SNR (dB)';
%ax.YLabel.String='Sum Rate (bits/sec/Hz)';
ax.YLabel.String='EE (bits/sec/Hz/W)';
ax.XGrid='on';
ax.YGrid='on';
ax.XMinorGrid = 'on';
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
lgd.FontSize = 8;
% 
% Figure Properties
ax.Units='centimeters';
ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
Fig.Units='centimeters';
Fig.Position = [10 5 (ax.TightInset(1)+ 2*ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];

%  Creat notation 
Not1=annotation('textbox');
X_points=[ .03*(diff(ax.XLim)) .3*(diff(ax.XLim))]+ax.XLim(1);
Y_points=[ .5*(diff(ax.YLim)) .67*(diff(ax.YLim))]+ax.YLim(1);
str1=['$N_t=$ ',num2str(Nt)];
str2=['$K=$',num2str(K)];
str3=['$N_{c}=$ ',num2str(Ncom)];
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