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
addpath('C:\Users\Jamal\OneDrive\Ph.D\Matlab\Functions') % Add path of the "Functions" folder.
addpath('C:\Users\Jamal\OneDrive\Ph.D\Matlab\Functions\Antenna Selection') % Add path of the "Functions" folder.

% Bits = 3;                      % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 8;                      % Number of the transmitter antennas on z axis
Txy = 8;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 8;                       % Number of the receiver antennas on y axis
Ns  = 6;                       % Number of data streams
NRF = 8;                      % Number of RF chains
Ncl =5;                       % Number of Channel clusters(Scatters)
Nray= 10;                       % Number of rays in each cluster
% AngSpread = 10/180*pi;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
SNR_dB = 0:5:25;             % Signal to noise ratio in dB
realization = 200;             % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_OMP=zeros(length(SNR_dB),realization);
R_Dig=zeros(length(SNR_dB),realization);
R_J=zeros(length(SNR_dB),realization);
R_Shrp=zeros(length(SNR_dB),realization);
R_DPS=zeros(length(SNR_dB),realization);
%%
F=64;
for r=1:length(SNR_dB)
    for reali=1:realization
        % generate channel
        [H,At,Ar,Alpha]=ChannelOFDM_MIMO(Tx,Rx,Ncl,'Nray',10,'F',F);
        % Full Digital beamforming
        [WOPT,FOPT]=FD_OFDM_MIMO(H,8,SNR(r));
        [WDp,WRFp,VRFp,VDp]=Sohrabi_Alg4_Partial(H,NRF,Ns,SNR(r));
        [WDs,WRFs,VRFs,VDs]=Sohrabi_Alg4_Full(H,NRF,Ns,SNR(r));
        [WD,WRF,VRF,VD]=PE_AltMin_OFDM_MIMO(H,NRF,Ns,SNR(r));
        for f=1:F 
            
            Fopt=FOPT(:,:,f);
            Wopt=WOPT(:,:,f);
%             trace(Fopt'*Fopt)
            
            Vtp=VRFp*VDp(:,:,f);           
            Wtp=WRFp*WDp(:,:,f);
%             trace(Vtp'*Vtp)
            
            Vts=VRFs*VDs(:,:,f);           
            Wts=WRFs*WDs(:,:,f);
%             trace(Vts'*Vts)
            
            Vt=VRF*VD(:,:,f);
            Wt=WRF*WD(:,:,f);
%             trace(Vt'*Vt)
           
            R_Dig(r,reali)  = R_Dig(r,reali) + log2(det(eye(8) +  pinv(Wopt) * H(:,:,f) * (Fopt) * Fopt'  * H(:,:,f)' * Wopt))/F;
            R_J(r,reali)  = R_J(r,reali) + log2(det(eye(Ns) +  pinv(Wt)   * H(:,:,f) *   Vt   * (Vt)'  * H(:,:,f)' * Wt))/F;
            R_Shrp(r,reali) = R_Shrp(r,reali)+ log2(det(eye(Ns) +  pinv(Wts)  * H(:,:,f) *   Vts  * (Vts)' * H(:,:,f)' * Wts))/F;
            R_DPS(r,reali) = R_DPS(r,reali)+ log2(det(eye(Ns) +  pinv(Wtp)  * H(:,:,f) *   Vtp  * (Vtp)' * H(:,:,f)' * Wtp))/F;

        end             
     end
end
% prepare data to plot
Xaxis=SNR_dB;

Data1= sum((R_Dig),2)/realization;
Data2= sum(R_J,2)/realization;
Data3= sum(R_DPS,2)/realization;
Data4= sum(R_Shrp,2)/realization;

%% plot 
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
P4=plot(Xaxis,Data4);

Leg1='Optimal Fully Digital BF';
Leg2='PE-AltMin';
Leg3='Alg. proposed by Sohrabi for Partially-connected';
Leg4='Alg. proposed by Sohrabi for Fully-connected';
% P1 properties
set(P1,'LineStyle','-');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','none');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color','k');
set(P1,'MarkerSize',7);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

% P2 properties
set(P2,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color',[0,.4,.4]);
set(P2,'MarkerSize',7);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% % P3 properties
set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','+');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color','b');
set(P3,'MarkerSize',7);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% % P3 properties
set(P4,'LineStyle','-.');           % '-' , '--' , ':' , '-.'
set(P4,'Marker','*');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P4,'LineWidth',2);
set(P4,'Color',[.7,0,0]);
set(P4,'MarkerSize',7);
P4_group= hggroup;set(P4,'Parent',P4_group);
set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');


lgd=legend(Leg1,Leg2,Leg3,Leg4);

% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
ax.YLim = [10 max(abs(Data1))];
ax.XTick = Xaxis;
ax.XTickLabel=Xaxis;
ax.XLabel.String='SNR (dB)';
ax.YLabel.String='Average Spectral Efficiency (bits/s/Hz)';
ax.XGrid='on';
ax.YGrid='on';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'on';
ax.TickLabelInterpreter='latex';
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

% % Figure Properties
% ax.Units='centimeters';
% ax.Position = [ax.TightInset(1) ax.TightInset(2) 11 8];
% Fig.Units='centimeters';
% Fig.Position = [10 5 (ax.TightInset(1)+ ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];




