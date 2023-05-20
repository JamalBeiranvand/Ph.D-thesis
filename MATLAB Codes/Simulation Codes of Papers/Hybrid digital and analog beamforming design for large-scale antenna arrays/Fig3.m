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
addpath('F:\Ph.D\Matlab\Functions') % Add path of the "Functions" folder.
b = 1;                      % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 1;                      % Number of the transmitter antennas on z axis
Txy = 10;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 10;                       % Number of the receiver antennas on y axis
Ns  = 2;                       % Number of data streams
NRF = Ns;                      % Number of RF chains
Ncl =15;                       % Number of Channel clusters(Scatters)
Nray= 1;                       % Number of rays in each cluster
% AngSpread = 10/180*pi;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
SNR_dB = 0:5:30;             % Signal to noise ratio in dB
realization = 100;             % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
Nq=2^b;
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_OMP=zeros(length(SNR_dB),realization);
R_Dig=zeros(length(SNR_dB),realization);
R_shrb=zeros(length(SNR_dB),realization);
R_shrbq=zeros(length(SNR_dB),realization);
%% Main code
for r=1:length(SNR_dB)
    for reali=1:realization
        % generate channel
        [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl);
        
        % Full Digital beamforming
        [W_Dig,F_Dig]=DigitalBeamforming(H,Ns);
        
        % Sihrabi Algorithm, quantize during computing
        [WDq,WRFq,VRFq,VDq]=Sohrabi_Alg(H,NRF,Ns,SNR(r),2^b);
        Vtq=VRFq*VDq;       
        Wtq=WRFq*WDq;
        
        % Sohrabi Algorithm, quantize after computing
        [WD,WRF,VRF,VD]=Sohrabi_Alg(H,NRF,Ns,SNR(r));
        VRF1=Qtz(VRF,Nq);
        WRF1=Qtz(WRF,Nq);
        Vt=VRF1*VD;
        Wt=WRF1*WD;
       
        
        % OMP Algorithm ,quantized
        [WRF_OPM,WBB_OPM,FRF_OMP, FBB_OMP] = OMP_Alg(H,At,Ar,NRF,Ns,SNR(r),Nq);

        
        R_shrbq(r,reali)  = log2(det(eye(Ns) +  pinv(Wt)  * H * Vt  * (Vt)'  * H' * Wt));
        R_shrb(r,reali) = log2(det(eye(Ns) +  pinv(Wtq)  * H * Vtq  * (Vtq)'  * H' * Wtq));        
        R_Dig(r,reali)   = log2(det(eye(Ns) + SNR(r)/Ns * pinv(W_Dig)  * H * F_Dig  * (F_Dig)'  * H' * W_Dig));
        R_OMP(r,reali)   = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_OPM * WBB_OPM)* H* FRF_OMP * (FBB_OMP) * FBB_OMP' * FRF_OMP' * H' * WRF_OPM * WBB_OPM));

    end
end
% prepare data to plot
Xaxis=SNR_dB;
Data1= sum((R_Dig),2)/realization;
Data2=sum(R_shrb,2)./realization;
Data3=sum(R_shrbq,2)./realization;
Data4= sum(R_OMP,2)/realization;
%% plot 
Fig=figure;
hold on
% P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
P4=plot(Xaxis,Data4);
%
% Legend1='Optimal Fully Digital BF';
Legend2='Proposed Alg. $b=1$';
Legend3='Quantized $-$ Proposed Alg. for $b=\infty$';
Legend4='Quantized $-$ Hybrid BF in [31] (OMP)';
% Legend5=' Quantized-Hybrid BF in [31] (OMP)';

% % P1 properties
% set(P1,'LineStyle',':');           % '-' , '--' , ':' , '-.' , 'none'
% set(P1,'Marker','none');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
% set(P1,'LineWidth',2);
% set(P1,'Color','k');
% set(P1,'MarkerSize',7);
% P1_group= hggroup;set(P1,'Parent',P1_group);
% set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% % P1.MarkerIndices =[ 1 6];

% P2 properties
set(P2,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color','b');
set(P2,'MarkerSize',9);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P3 properties
set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','s');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color','b');
set(P3,'MarkerSize',10);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P3 properties
set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P4,'Marker','v');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P4,'LineWidth',2);
set(P4,'Color','g');
set(P4,'MarkerSize',7);
P4_group= hggroup;set(P4,'Parent',P4_group);
set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

lgd=legend(Legend2,Legend3,Legend4);

% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
ax.YLim = [4 26];
ax.XTick = Xaxis;
ax.XTickLabel=Xaxis;
ax.YTick = 4:2:26;
ax.XLabel.String='SNR (dB)';
ax.YLabel.String='Spectral Efficiency (bits/sec/Hz)';
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

