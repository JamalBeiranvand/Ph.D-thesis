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
Bits = 'inf';                      % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz =[8,16];                      % Number of the transmitter antennas on z axis
Txy =[8,16];                      % Number of the transmitter antennas on y axis
Rxz =[4,8];                       % Number of the receiver antennas on z axis
Rxy =[4,8];                       % Number of the receiver antennas on y axis
Ns  = 1;                       % Number of data streams
NRF = 4;                      % Number of RF chains
Ncl =8;                       % Number of Channel clusters(Scatters)
Nray= 10;                       % Number of rays in each cluster
AngSpread =0:15;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
TxSector=[60 360];
RxSector=[20 360];
SNR_dB = 0;             % Signal to noise ratio in dB
realization = 500;             % Iteration of simulation
%% Initialize Parameters
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_opt=zeros(length(SNR_dB),realization);    
R_OMP=zeros(length(SNR_dB),realization);

Data1=zeros(length(Txz),length(AngSpread));
Data2=zeros(length(Txz),length(AngSpread));
%% Main code 
for t=1:length(Txz)
    Tx  = [Txz(t),Txy(t)];               % Srtucture of the transmit antanna array
    Rx  = [Rxz(t),Rxy(t)];               % Srtucture of the receiver antanna array 
    Nt  = Txz(t)*Txy(t);                 % Number of the transmit antennas
    Nr  = Rxz(t)*Rxy(t);                 % Number of the receive antennas

    for n = 1:length(AngSpread)
        for reali = 1:realization
            % generate channel
             [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl,'Nray',Nray,'AngSpread',AngSpread(n),'TxSector',TxSector,'RxSector',RxSector );

             for r=1:length(SNR)
                 %Full Digital beamforming
                [Wopt,Fopt]=DigitalBeamforming(H,Ns);

                % obtain FRF and FBB by function of OMP algorithm
                [WRF_OPM,WBB_OPM,FRF_OMP, FBB_OMP] = OMP_Alg(H,At,Ar,NRF,Ns,SNR(r));

                 % compute Spectral efficiency for each algorithm 
                % compute Optimal Spectral efficiency         
                 R_opt(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(Wopt) * H * Fopt * (Fopt)' * H' * Wopt));
                 R_OMP(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_OPM * WBB_OPM)* H * FRF_OMP * (FBB_OMP) * FBB_OMP' * FRF_OMP' * H' * WRF_OPM * WBB_OPM));
             end
        end  
        Data1(t,n)= sum(R_opt,2)/realization;
        Data2(t,n)= sum(R_OMP,2)/realization;
    end
end

Tx  = [8,8];             % Srtucture of the transmit antanna array
Rx  = [4,4];               % Srtucture of the receiver antanna array 
Ns=2;
for n = 1:length(AngSpread)
    for reali = 1:realization
        % generate channel
         [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl,'Nray',Nray,'AngSpread',AngSpread(n),'TxSector',TxSector,'RxSector',RxSector );

         for r=1:length(SNR)
             %Full Digital beamforming
            [Wopt,Fopt]=DigitalBeamforming(H,Ns);

            % obtain FRF and FBB by function of OMP algorithm
            [WRF_OPM,WBB_OPM,FRF_OMP, FBB_OMP] = OMP_Alg(H,At,Ar,NRF,Ns,SNR(r));

             % compute Spectral efficiency for each algorithm 
            % compute Optimal Spectral efficiency         
             R_opt(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(Wopt) * H * Fopt * (Fopt)' * H' * Wopt));
             R_OMP(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_OPM * WBB_OPM)* H * FRF_OMP * (FBB_OMP) * FBB_OMP' * FRF_OMP' * H' * WRF_OPM * WBB_OPM));
         end
    end  
    Data1(3,n)= sum(R_opt,2)/realization;
    Data2(3,n)= sum(R_OMP,2)/realization;
end

NRF=6;
for n = 1:length(AngSpread)
    for reali = 1:realization
        % generate channel
         [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl,'Nray',Nray,'AngSpread',AngSpread(n),'TxSector',TxSector,'RxSector',RxSector );

         for r=1:length(SNR)
            % obtain FRF and FBB by function of OMP algorithm
            [WRF_OPM,WBB_OPM,FRF_OMP, FBB_OMP] = OMP_Alg(H,At,Ar,NRF,Ns,SNR(r));

             % compute Spectral efficiency for each algorithm         
             R_OMP(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_OPM * WBB_OPM)* H * FRF_OMP * (FBB_OMP) * FBB_OMP' * FRF_OMP' * H' * WRF_OPM * WBB_OPM));
         end
    end  
    Data2(4,n)= sum(R_OMP,2)/realization;
end
% prepare Data
Xaxis=AngSpread;
%% plot 
Fig=figure;
hold on
P1=plot(Xaxis,Data1);
P2=plot(Xaxis,Data2);
% P3=plot(Xaxis,Data3);
% P4=plot(Xaxis,Data4);
% P5=plot(Xaxis,Data5);

Legend1='Optimal Digital Beamforming';
Legend2='Sparse Precoding & Combining (OMP Algorithm) ';
% Legend3='OMP Algorithm';
% Legend4='DPS structure';
% Legend5='PE-AltMin Algorithm';

% P1 properties
set(P1,'LineStyle','-');           % '-' , '--' , ':' , '-.' , 'none'
set(P1,'Marker','S');       % 'o' , '+'  , '*' , 'x' ,'s' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram', 'none'
set(P1,'LineWidth',2);
set(P1,'Color',[0 .7 .2]);
set(P1,'MarkerSize',10);
P1_group= hggroup;set(P1,'Parent',P1_group);
set(get(get(P1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
% P1.MarkerIndices =[ 1 6];

% P2 properties
set(P2,'LineStyle','-');           % '-' , '--' , ':' , '-.'
set(P2,'Marker','pentagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
set(P2,'LineWidth',2);
set(P2,'Color',[0 .2 .8]);
set(P2,'MarkerSize',5);
P2_group= hggroup;set(P2,'Parent',P2_group);
set(get(get(P2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% % P2 properties
% set(P3,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P3,'Marker','hexagram');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P3,'LineWidth',2);
% set(P3,'Color',[.6 .2 .2]);
% set(P3,'MarkerSize',7);
% P3_group= hggroup;set(P3,'Parent',P3_group);
% set(get(get(P3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% % 
% % P2 properties
% set(P4,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P4,'Marker','^');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P4,'LineWidth',2);
% set(P4,'Color',[0 .6 .2]);
% set(P4,'MarkerSize',7);
% P4_group= hggroup;set(P4,'Parent',P4_group);
% set(get(get(P4_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% 
% % P2 properties
% set(P5,'LineStyle','-');           % '-' , '--' , ':' , '-.'
% set(P5,'Marker','x');       % 'o' , '+'  , '*' , 'x' , 's' , '^' , 'v' , '>' , '<' , 'pentagram', 'hexagram'
% set(P5,'LineWidth',2);
% set(P5,'Color',[.5 .5 0]);
% set(P5,'MarkerSize',7);
% P5_group= hggroup;set(P5,'Parent',P5_group);
% set(get(get(P5_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');


lgd=legend(Legend1,Legend2);
% Axes Properties 
ax=gca;
ax.XLim = [min(Xaxis) max(Xaxis)];
 ax.YLim = [8 15];
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

