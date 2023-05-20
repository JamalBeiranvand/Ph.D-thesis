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
Txy = 64;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 16;                       % Number of the receiver antennas on y axis
Ns  = 4;                       % Number of data streams
NRF = 5;                      % Number of RF chains
Ncl =15;                       % Number of Channel clusters(Scatters)
Nray= 1;                       % Number of rays in each cluster
% AngSpread = 10/180*pi;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
SNR_dB = -10:2:6;             % Signal to noise ratio in dB
realization = 200;             % Iteration of simulation
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
R_shrbq1=zeros(length(SNR_dB),realization);
R_shrbq3=zeros(length(SNR_dB),realization);
R_shrbq=zeros(length(SNR_dB),realization);
%% Main code
for r=1:length(SNR_dB)
    for reali=1:realization
        % generate channel
        [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl);
        
        % Full Digital beamforming
        [W_Dig,F_Dig]=DigitalBeamforming(H,Ns);
        
        % Sohrabi Algorithm,
        [WD,WRF,VRF,VD]=Sohrabi_Alg(H,NRF,Ns,SNR(r));        
        Vt=VRF*VD;
        Wt=WRF*WD;
        
        % Sohrabi Algorithm, NRF =NS
        [WDq,WRFq,VRFq,VDq]=Sohrabi_Alg(H,Ns,Ns,SNR(r),Nq);
        Vtq=VRFq*VDq;       
        Wtq=WRFq*WDq;
        
        % Sohrabi Algorithm, NRF =NS+1
        [WDq1,WRFq1,VRFq1,VDq1]=Sohrabi_Alg(H,Ns+1,Ns,SNR(r),Nq);
        Vtq1=VRFq1*VDq1;       
        Wtq1=WRFq1*WDq1;
        
        % Sohrabi Algorithm, NRF =NS+3
        [WDq3,WRFq3,VRFq3,VDq3]=Sohrabi_Alg(H,Ns+3,Ns,SNR(r),Nq);
        Vtq3=VRFq3*VDq3;       
        Wtq3=WRFq3*WDq3;
        
        % Sohrabi Algorithm, quantize after computing
        [WD,WRF,VRF,VD]=Sohrabi_Alg(H,NRF,Ns,SNR(r));
        VRF1=Qtz(VRF,Nq);
        WRF1=Qtz(WRF,Nq);
        VtQ=VRF1*VD;
        WtQ=WRF1*WD;
        
        R_Dig(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(W_Dig) * H * F_Dig * (F_Dig)' * H' * W_Dig));
        R_shrb(r,reali) = log2(det(eye(Nr) + Wt*inv((Wt)'*Wt)*(Wt)' * H * Vt * (Vt)' * H'));
        R_shrbq(r,reali) = log2(det(eye(Nr) + Wtq*inv((Wtq)'*Wtq)*(Wtq)' * H * Vtq * (Vtq)' * H'));
        R_shrbq1(r,reali) = log2(det(eye(Nr) + Wtq1*inv((Wtq1)'*Wtq1)*(Wtq1)' * H * Vtq1 * (Vtq1)' * H'));
        R_shrbq3(r,reali) = log2(det(eye(Nr) + Wtq3*inv((Wtq3)'*Wtq3)*(Wtq3)' * H * Vtq3 * (Vtq3)' * H'));
        R_shrbQ(r,reali) = log2(det(eye(Nr) + WtQ*inv((WtQ)'*WtQ)*(WtQ)' * H * VtQ * (VtQ)' * H'));


    end
end
R_shrbq(isnan(R_shrbq))=0;
R_shrbq1(isnan(R_shrbq1))=0;
R_shrbq3(isnan(R_shrbq3))=0;
R_shrbQ(isnan(R_shrbQ))=0;
R_shrb(isnan(R_shrb))=0;
%% 
Data1= sum((R_Dig),2)/realization;
Data2= sum(R_shrb,2)   ./(realization-sum(R_shrb==0,2));
Data3= sum(R_shrbq,2)  ./(realization-sum(R_shrbq==0,2));
Data4= sum(R_shrbq1,2) ./(realization-sum(R_shrbq1==0,2));
Data5= sum(R_shrbq3,2) ./(realization-sum(R_shrbq3==0,2));
Data6= sum(R_shrbQ,2) ./(realization-sum(R_shrbQ==0,2));
Xaxis= SNR_dB;
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
P3=plot(Xaxis,Data3);
P4=plot(Xaxis,Data4);
P5=plot(Xaxis,Data5);
P6=plot(Xaxis,Data6);

Legend1='Optimal Fully Digital BF';
Legend2='Proposed Alg for $b=\infty$, $N_{RF}=N_s $';
Legend3='Proposed Alg for $b=1$, $N_{RF}=N_s $';
Legend4='Proposed Alg for $b=1$, $N_{RF}=N_s+1 $';
Legend5='Proposed Alg for $b=1$, $N_{RF}=N_s+3 $';
Legend6='Quantized$-$Proposed  Alg for $b=\infty$';
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
set(P2,'LineStyle','-.');           % '-' , '--' , ':' , '-.'
set(P2,'Marker', 'o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color','k');
set(P2,'MarkerSize',9);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% P2 properties
set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P3,'Marker','o');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P3,'LineWidth',2);
set(P3,'Color',[0 0 .5]);
set(P3,'MarkerSize',7);
P3_group= hggroup;set(P3,'Parent',P3_group);
set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% P2 properties
set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P4,'Marker','s');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P4,'LineWidth',2);
set(P4,'Color',[0 0 .7]);
set(P4,'MarkerSize',9);
P4_group= hggroup;set(P4,'Parent',P4_group);
set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P2 properties
set(P5,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P5,'Marker','>');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P5,'LineWidth',2);
set(P5,'Color',[0 0 .9]);
set(P5,'MarkerSize',7);
P5_group= hggroup;set(P5,'Parent',P5_group);
set(get(get(P5_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

% P2 properties
set(P6,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P6,'Marker','v');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P6,'LineWidth',2);
set(P6,'Color',[.7 0 0]);
set(P6,'MarkerSize',7);
P6_group= hggroup;set(P6,'Parent',P6_group);
set(get(get(P6_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');


lgd=legend(Legend1,Legend2,Legend3,Legend4,Legend5,Legend6);
% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
ax.YLim = [4 36];
ax.XTick = Xaxis;
ax.YTick = 4:4:36;
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
% % % % 

