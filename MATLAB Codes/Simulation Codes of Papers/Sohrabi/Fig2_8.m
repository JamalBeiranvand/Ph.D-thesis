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
addpath(genpath('F:\Ph.D\Matlab\Functions') )% Add path of the "Functions" folder.
K=8;                           % number of users
% Bits = 3;                    % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 1;                       % Number of the transmitter antennas on z axis
Txy = 64;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 1;                       % Number of the receiver antennas on y axis
% Ns  = 1;                       % Number of data streams
NrfTx = 9;                     % Number of RF chains at BS
Ncl =15;                       % Number of Channel clusters(Scatters)
Nray= 1;                       % Number of rays in each cluster
SNR_dB = -10:4:10;             % Signal to noise ratio in dB
realization = 20;              % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_ZF=zeros(length(SNR_dB),realization);
R_Shr=zeros(length(SNR_dB),realization);
R=zeros(K,1);
%%
for r=1:length(SNR_dB)
    disp(r)
    Snr=SNR(r);
    for reali=1:realization
        % generate channel
        [H,At,Ar,~]=MU_Channel(Tx,K,Ncl); 
        
        % Full Digital ZF beamforming
        AZF=ZF_Down(H);
        G= diag(ones(K,1))*sqrt(Snr);
        R_ZF(r,reali)=SumRate(H,AZF,G);
        
        % Sohrabi Algorithm
        [VRF,VD]= Sohrabi_Alg3(H,NrfTx,Snr);
        sum(sum(abs(VRF*VD).^2))
        R_Shr(r,reali)=SumRate(H,VRF*VD,eye(K));
    end
end
Data1= sum((R_ZF),2)/realization;
Data2= sum(R_Shr,2)/realization;
% Data3= sum(R_OMP,2)/realization;
%% plot 

Xaxis= SNR_dB;
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
% P3=plot(Xaxis,Data3);
% P4=plot(Xaxis,Data4);
% P5=plot(Xaxis,Data5);

Legend1='Fully Digital ZF';
Legend2='Proposed Algorithm $N^{RF}=9$';
% Legend3='Proposed Hybrid BF Alg for $b=1$, $N_{RF}=N_s+1 $';
% Legend4='Proposed Hybrid BF Alg for $b=1$, $N_{RF}=N_s+3 $';
% Legend5=' Quantized-Hybrid BF in [31] (OMP)';

% P1 properties
set(P1,'LineStyle',':');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','none');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color','k');
set(P1,'MarkerSize',10);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

% P2 properties
set(P2,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color',[0 0 1]);
set(P2,'MarkerSize',5);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% % P2 properties
% set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P3,'Marker','hexagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P3,'LineWidth',2);
% set(P3,'Color',[0 0 .6]);
% set(P3,'MarkerSize',7);
% P3_group= hggroup;set(P3,'Parent',P3_group);
% set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% % 
% % P2 properties
% set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P4,'Marker','^');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P4,'LineWidth',2);
% set(P4,'Color',[.3 0 0]);
% set(P4,'MarkerSize',7);
% P4_group= hggroup;set(P4,'Parent',P4_group);
% set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% % P2 properties
% set(P5,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P5,'Marker','x');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P5,'LineWidth',2);
% set(P5,'Color',[0 .5 0]);
% set(P5,'MarkerSize',7);
% P5_group= hggroup;set(P5,'Parent',P5_group);
% set(get(get(P5_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');


lgd=legend(Legend1,Legend2);
% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
ax.YLim = [0 50];
ax.XTick = Xaxis;
ax.YTick = 0:5:50;
ax.XTickLabel=Xaxis;
ax.XLabel.String='SNR (dB)';
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

