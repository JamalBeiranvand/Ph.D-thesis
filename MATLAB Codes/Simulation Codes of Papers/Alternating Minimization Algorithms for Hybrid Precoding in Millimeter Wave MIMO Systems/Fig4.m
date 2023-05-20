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
% Bits = 3;                      % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 12;                      % Number of the transmitter antennas on z axis
Txy = 12;                      % Number of the transmitter antennas on y axis
Rxz = 6;                       % Number of the receiver antennas on z axis
Rxy = 6;                       % Number of the receiver antennas on y axis
Ns  = 2;                       % Number of data streams
NRF_Vec = 2:6;                      % Number of RF chains
Ncl =5;                       % Number of Channel clusters(Scatters)
Nray= 10;                       % Number of rays in each cluster
AngSpread = 10;              % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
SNR_dB = 0;             % Signal to noise ratio in dB
realization = 20;             % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_MO=zeros(length(SNR_dB),realization);    
R_Dig=zeros(length(SNR_dB),realization);
R_SDR=zeros(length(SNR_dB),realization);
R_OMP=zeros(length(SNR_dB),realization);

% Data1=zeros(length(SNR_dB),length(NRF_Vec));
% Data2=zeros(length(SNR_dB),length(NRF_Vec));
% Data3=zeros(length(SNR_dB),length(NRF_Vec));
% Data4=zeros(length(SNR_dB),length(NRF_Vec));
%% Main code
% for n=1:length(NRF_Vec)
%     NRF = NRF_Vec(n);
%     for reali=1:realization
%         disp(((n-1)*reali+reali)/(realization*length(NRF_Vec))*100)
%         % generate channel
%         [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl,'Nray',Nray,'AngSpread',AngSpread);
% 
%         % Full Digital beamforming
%         [W_opt,F_opt]=DigitalBeamforming(H,Ns);
% 
%         % MO-AltMin Algorithm
%         [WRF_MO,WBB_MO,FRF_MO,FBB_MO] = MO_AltMin(H,NRF,Ns);
%         
%         if (NRF~=5)
%            % SDR-AltMin Algorithm
%            [FRF_SDR,FBB_SDR] = PartSDR_AltMin(H,NRF,Ns);
%         end
%         
%        %OMP algorithm
%        [WRF_OPM,WBB_OPM,FRF_OMP,FBB_OMP] = OMP_Alg(H,At,Ar,NRF,Ns,SNR);
% 
%         R_OMP(n,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(WRF_OPM * WBB_OPM)  * H * FRF_OMP * FBB_OMP * (FBB_OMP)' * FRF_OMP' * H' * WRF_OPM * WBB_OPM));
%         R_SDR(n,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(W_opt)  * H * FRF_SDR * FBB_SDR * (FBB_SDR)' * FRF_SDR' * H' * W_opt));
% 
%         R_MO(n,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(WRF_MO * WBB_MO)  * H * FRF_MO * FBB_MO * (FBB_MO)' * FRF_MO' * H' * WRF_MO * WBB_MO));
%         R_Dig(n,reali) = log2(det(eye(Ns) + SNR/Ns * pinv(W_opt)  * H * F_opt  * (F_opt)'  * H' * W_opt));
%     end
% end
% Data1= sum(R_Dig,2)/realization;
% Data1(:)=mean(Data1);
% Data2= sum(R_MO,2)/realization;
% Data3= sum(R_OMP,2)/realization;
% Data4= sum(R_SDR,2)/realization;

% % prepare data to plot
% DataFig4(:,1)=Data1;
% DataFig4(:,2)=Data2;
% DataFig4(:,3)=Data3;
% DataFig4(:,4)=Data4;
% save('DataFig4')

Xaxis=NRF_Vec;
%% plot 

load('DataFig4')
Data1=DataFig4(:,1);
Data2=DataFig4(:,2);
Data3=DataFig4(:,3);
Data4=DataFig4(:,4);


Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
P4=plot([2:4,6],Data4([1:3,5],1));
% P5=plot(Xaxis,Data5);
% P6=plot(Xaxis,Data6);

Legend1='Optimal Digital Beamforming';
Legend2='MO-AltMin';
Legend3='OMP Algorithm';
Legend4='SDR-AltMin';
% Legend5='SIC-Based Method';
% Legend6='Analog Beamforming';

% P1 properties
set(P1,'LineStyle','-');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','o');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color','r');
set(P1,'MarkerSize',10);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

% P2 properties
set(P2,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','pentagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color',[.7 0 1]);
set(P2,'MarkerSize',5);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P2 properties
set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','^');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color',[0 .5 0]);
set(P3,'MarkerSize',7);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% P2 properties
set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P4,'Marker','diamond');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P4,'LineWidth',2);
set(P4,'Color',[0.8 0.5 0]);
set(P4,'MarkerSize',7);
P4_group= hggroup;set(P4,'Parent',P4_group);
set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% % P2 properties
% set(P5,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P5,'Marker','+');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P5,'LineWidth',2);
% set(P5,'Color',[.1 .2 .5]);
% set(P5,'MarkerSize',10);
% P5_group= hggroup;set(P5,'Parent',P5_group);
% set(get(get(P5_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% % P6 properties
% set(P6,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P6,'Marker','*');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P6,'LineWidth',2);
% set(P6,'Color',[.3 .2 .5]);
% set(P6,'MarkerSize',10);
% P6_group= hggroup;set(P6,'Parent',P6_group);
% set(get(get(P6_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 


lgd=legend(Legend1,Legend2,Legend3,Legend4);
% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
ax.YLim = [13.5 19];
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
% % Figure Properties
% ax.Units='centimeters';
% ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
% Fig.Units='centimeters';
% Fig.Position = [10 5 (ax.TightInset(1)+ 2*ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];

% %  Creat notation 
% Not1=annotation('textbox');
% X_points=[-28 0];
% Y_points=[ .37*(diff(ax.YLim)) .65*(diff(ax.YLim))];
% str1=['UPA structure at Tx: ',num2str(Tx(1)),'$\times$',num2str(Tx(2))];
% str2=['UPA structure at Rx: ',num2str(Rx(1)),'$\times$',num2str(Rx(2))];
% str3=['$N_{s}=$ ',num2str(Ns),',  $N_{RF}=$ ',num2str(NRF)];
% str4=['$N_{cl}=$ ',num2str(Ncl),',  $N_{ray}=$ ',num2str(Nray)];
% str5=['$b=$',num2str(b)];
% Not1.String={str1,str2,str3,str4, str5};
% Not1.FontSize = 10;
% Not1.Interpreter='latex';
% Not1.Units='centimeters';
% Dx=ax.Position(3)/(diff(ax.XLim));
% Dy=ax.Position(4)/(diff(ax.YLim));
% x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% Not1.Position=[x_begin y_begin  x_End y_End];
% % % % % 
