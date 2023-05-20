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
addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' ))% Add path of the "Functions" folder.

Bits = 3;                      % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 12;                      % Number of the transmitter antennas on z axis
Txy = 12;                      % Number of the transmitter antennas on y axis
Rxz = 6;                       % Number of the receiver antennas on z axis
Rxy = 6;                       % Number of the receiver antennas on y axis
Ns  = 20;                       % Number of data streams
NRF = Ns;                      % Number of RF chains
Ncl =5;                       % Number of Channel clusters(Scatters)
Nray= 3;                       % Number of rays in each cluster
AngSpread = 10;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
SNR_dB = -30:5:10;             % Signal to noise ratio in dB
realization = 100;             % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_Nov=zeros(length(SNR_dB),realization);  
R_NovP=zeros(length(SNR_dB),realization); 
R_DigP=zeros(length(SNR_dB),realization);

Data1=zeros(length(SNR_dB),length(Bits));
Data2=zeros(length(SNR_dB),length(Bits));
Data3=zeros(length(SNR_dB),length(Bits));


%% Main code
for n=1:length(Bits)   
    b=Bits(n);
    Nps=2^b;
    % Generate feasable Set
    [Setb,Sb]=GenUniqueSet(Nps);
    
    for reali=1:realization
        
        % generate channel
        [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl,'Nray',Nray,'AngSpread',AngSpread);
        
        % Full Digital beamforming
        [W_Dig,F_Dig,lambda]=DigitalBeamforming(H,Ns);
        
        % Novel Algorithm
        [WRF_Nov,WBB_Nov,FRF_Nov,FBB_Nov] = NovelAlg(H,Ns,Setb,Sb);
          

          for r=1:length(SNR)
%               g=WaterFilling(lambda,Nt,SNR(r)); 
              g=WaterFilling1(lambda,SNR(r)); 
              G=diag(g);
              R_Nov(r,reali)  = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_Nov * WBB_Nov)  * H * FRF_Nov * FBB_Nov * (FBB_Nov)' * FRF_Nov' * H' * WRF_Nov * WBB_Nov));
              R_NovP(r,reali) = log2(det(eye(Ns) +      G    * pinv(WRF_Nov * WBB_Nov)  * H * FRF_Nov * FBB_Nov * (FBB_Nov)' * FRF_Nov' * H' * WRF_Nov * WBB_Nov));
              R_DigP(r,reali) = log2(det(eye(Ns) +      G    * pinv(W_Dig)  * H * F_Dig  * (F_Dig)'  * H' * W_Dig));
          end
    end
    Data1(:,n)= sum(R_DigP,2)/realization;
    Data2(:,n)= sum(R_NovP,2)/realization;
    Data3(:,n)= sum(R_Nov,2)/realization;
end
% prepare data to plot
Xaxis=SNR_dB;
%% plot 
Fig=figure;
hold on
P1=plot(Xaxis,Data1(:,1));
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
% P4=plot(Xaxis,Data4);
% P5=plot(Xaxis,Data5);

Legend1='Optimal Digital Beamforming with power allocation';
Legend2=['Novel algorithm with power allocation for ','$b=$',num2str(b)];
Legend3=['Novel algoritthm without power allocation for ','$b=$',num2str(b)];
% Legend4='DPS structure';
% Legend5='PE-AltMin Algorithm';
% Legend3=('Data1 $N_t=144$, $Nr=36$');

% P1 properties
set(P1,'LineStyle','-');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','none');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color','r');
set(P1,'MarkerSize',10);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

% P2 properties
set(P2,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','hexagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color',[0 .6 .5]);
set(P2,'MarkerSize',12);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% % P2 properties
set(P3,'LineStyle','-.');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','pentagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color',[0 .2 .6]);
set(P3,'MarkerSize',10);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% % 
% % P2 properties
% set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P4,'Marker','hexagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P4,'LineWidth',2);
% set(P4,'Color',[.8 .5 0]);
% set(P4,'MarkerSize',7);
% P4_group= hggroup;set(P4,'Parent',P4_group);
% set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% % P2 properties
% set(P5,'LineStyle','--');           % '-' , '--' , ':' , '-.'
% set(P5,'Marker','x');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P5,'LineWidth',2);
% set(P5,'Color',[.8 .5 .5]);
% set(P5,'MarkerSize',7);
% P5_group= hggroup;set(P5,'Parent',P5_group);
% set(get(get(P5_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 

lgd=legend(Legend1,Legend2,Legend3);
% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
% ax.YLim = [5 40];
ax.XTick = Xaxis;
ax.XTickLabel=Xaxis;
ax.XLabel.String='SNR (dB)';
ax.YLabel.String='Spectral Efficiency (bits/sec/Hz)';
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

