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

K=10;                           % number of users
NcomM=K:30:100;
Nps=11;
Txz =16;                       % Number of the transmitter antennas on z axis
Txy =16 ;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 1;                       % Number of the receiver antennas on y axis
% Ns  = 1;                       % Number of data streams
% NrfTx = 9;                     % Number of RF chains at BS
Ncl =10;                       % Number of Channel clusters(Scatters)
Nray= 5;                       % Number of rays in each cluster
AngSpread = 10;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
TxSector=[360 360];
RxSector=[361 360];
SNR_dB = 10;             % Signal to noise ratio in dB
realization = 2000;              % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);  
%% Allocate memory space for matrices
R_Implemented=zeros(length(NcomM),realization);
R_ZF=zeros(length(NcomM),realization);
R_MRC=zeros(length(NcomM),realization);
SumR_Dynamic=zeros(length(NcomM),realization);
SumR_Fixed=zeros(length(NcomM),realization);
%%
[Setb,Sb]=GenUniqueSet(Nps);
% Fig=figure;
Rcell=cell(length(NcomM),1);

for r=1:length(NcomM)
    disp(r)
    Ncom=NcomM(r);
    Nu=Nt-Ncom;
    Snr=SNR;
    for reali=1:realization
        % generate channel
        [H,At,Ar,~]=MU_Channel(Tx,K,Ncl,'Nray',Nray,'AngSpread',AngSpread);
%         H=sqrt(1/2)*(randn(K,Nt)+1j*randn(K,Nt));

        G= diag(ones(K,1))*sqrt(Snr);        
        % Full Digital ZF beamforming
        Fzf=ZF_Down(H);        
        R_ZF(r,reali)=SumRate(H,Fzf,G);
        
        % Full Digital MRC beamforming
        Fmrc=MRC_Down(H);
        R_MRC(r,reali)=SumRate(H,Fmrc,G);
        
        H=H.';
         %% Fixed
        [KSets_Fixed,~]=NonSharedFixed(Nu,K);
        F_Fixed=SPF(H,KSets_Fixed);       
        [S,C]=BeamSwitch(F_Fixed,Setb,Sb);
        FRF_FSP=S*C;
        FRF_FSP=FRF_FSP/sqrt(sum(sum(abs(FRF_FSP).^2)));
        SumR_Fixed(r,reali)=SumRate(H.',FRF_FSP,G);
    end
end

%% plot 
Fig=figure;

for i=1:length(NcomM)
    subplot(length(NcomM),1,i)
    hist([R_ZF(i,:).',SumR_Fixed(i,:).'],200)
    ax=gca;
    ax.XLim = [0  max(max(R_ZF))+1];
    ax.XLabel.String=['Sum Rate (bits/Sec/Hz) for $N_{com}=$',num2str(NcomM(i))];
    ax.XLabel.Interpreter='latex';  
    ax.YLabel.String=['Counts'];
    ax.YLabel.Interpreter='latex'; 
end
Legend1='Fully-connected Architecture []';
%Legend2='Fully Digital MRC Precoder';
% Legend2='Dynamic Semi-Partial Architecture';
Legend3='Fixed Semi-Partial Architecture';
lgd=legend(Legend1,Legend3);


% Legend4='Full Search';
% Legend1='Fully Digital ZF Precoder';
% Legend2=['Dynamic Structure'];
% Legend3=['Fixed Structure '];
% 
% 
% lgd=legend(Legend1,Legend2,Legend3,Legend4);
% Xaxis= NcomM;
% figure(Fig);
% ax=gca;
% ax.XLim = [min(Xaxis)-1 max(Xaxis)+1];
% ax.YLim = [0 max(R_ZF)+20];
% ax.XTick = Xaxis;
% % ax.YTick = 0:5:50;
% ax.XTickLabel=Xaxis;
% ax.XLabel.String='$N_{com}$';
% ax.YLabel.String='Sum Rate (bits/sec/Hz)';
% ax.XGrid='on';
% ax.YGrid='on';
% ax.XMinorGrid = 'off';
% ax.YMinorGrid = 'on';
% ax.YLabel.Interpreter='latex';
% ax.XLabel.Interpreter='latex';
% ax.FontSize = 11;
% ax.LineWidth = .9;
% ax.Box='on';
% 
% % Legend Properties %
lgd.Location= 'northwest'; %%
lgd.Box='on';
lgd.LineWidth=.9;
lgd.Interpreter='latex';
lgd.FontSize = 10;
% % 
% % Figure Properties
% ax.Units='centimeters';
% ax.Position = [ax.TightInset(1)+ax.TightInset(3) ax.TightInset(2) 11 8];
% Fig.Units='centimeters';
% Fig.Position = [10 5 (ax.TightInset(1)+ 2*ax.TightInset(3)+ ax.Position(3) ) (ax.TightInset(2)+ ax.TightInset(4)+ ax.Position(4))];
% % 
% % %%  Creat notation 
% Not1=annotation('textbox');
% X_points=[0.6*diff(ax.XLim) .85*(diff(ax.XLim))]+ax.XLim(1);
% Y_points=[.73*(diff(ax.YLim)) .93*(diff(ax.YLim))]+ax.YLim(1);
% str1=['$N_t=$ ',num2str(Nt)];
% str2=['$K=$',num2str(K)];
% str3=['$SNR=$ ',num2str(SNR_dB),' dB'];
% Not1.String={str1,str2,str3};
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
% % Ell=annotation('ellipse');
% % X_points=[19 21];
% % Y_points=[25 55];
% % Ell.Units='centimeters';
% % Dx=ax.Position(3)/(diff(ax.XLim));
% % Dy=ax.Position(4)/(diff(ax.YLim));
% % x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% % y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% % x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% % y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% % Ell.Position=[x_begin y_begin  x_End y_End];
% % 
% % Textarrow1=annotation('textarrow');
% % X_points=[8 8];
% % Y_points=[30 62];
% % str1=['$b_{ps}^{t}=b_{ps}^{r}=\infty$'];
% % Textarrow1.String={str1};
% % Textarrow1.Interpreter='latex';
% % Textarrow1.FontSize = 10;
% % Textarrow1.Units='centimeters';
% % Dx=ax.Position(3)/(diff(ax.XLim));
% % Dy=ax.Position(4)/(diff(ax.YLim));
% % x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% % y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% % x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% % y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% % Textarrow1.Position=[x_begin y_begin  x_End y_End];
% % %%
% % Ell2=annotation('ellipse');
% % X_points=[19 21];
% % Y_points=[90 125];
% % Ell2.Units='centimeters';
% % Dx=ax.Position(3)/(diff(ax.XLim));
% % Dy=ax.Position(4)/(diff(ax.YLim));
% % x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% % y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% % x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% % y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% % Ell2.Position=[x_begin y_begin  x_End y_End];
% % 
% % Textarrow2=annotation('textarrow');
% % X_points=[16 16];
% % Y_points=[30 81];
% % str1=['$b_{ps}^{t}=b_{ps}^{r}=$',num2str(B_Vec{2})];
% % Textarrow2.String={str1};
% % Textarrow2.Interpreter='latex';
% % Textarrow2.FontSize = 10;
% % Textarrow2.Units='centimeters';
% % Dx=ax.Position(3)/(diff(ax.XLim));
% % Dy=ax.Position(4)/(diff(ax.YLim));
% % x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% % y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% % x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% % y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% % Textarrow2.Position=[x_begin y_begin  x_End y_End];
% % %%
% % Ell3=annotation('ellipse');
% % X_points=[19 21];
% % Y_points=[140 185];
% % Ell3.Units='centimeters';
% % Dx=ax.Position(3)/(diff(ax.XLim));
% % Dy=ax.Position(4)/(diff(ax.YLim));
% % x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% % y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% % x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% % y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% % Ell3.Position=[x_begin y_begin  x_End y_End];
% % 
% % Textarrow3=annotation('textarrow');
% % X_points=[24 24];
% % Y_points=[30 90];
% % str1=['$b_{ps}^{t}=b_{ps}^{r}=$',num2str(B_Vec{3})];
% % Textarrow3.String={str1};
% % Textarrow3.Interpreter='latex';
% % Textarrow3.FontSize = 10;
% % Textarrow3.Units='centimeters';
% % Dx=ax.Position(3)/(diff(ax.XLim));
% % Dy=ax.Position(4)/(diff(ax.YLim));
% % x_begin=ax.Position(1)-ax.XLim(1)*Dx+X_points(1)*Dx;
% % y_begin=ax.Position(2)-ax.YLim(1)*Dy+Y_points(1)*Dy;
% % x_End=(ax.Position(1)-ax.XLim(1)*Dx+X_points(2)*Dx)-x_begin;
% % y_End=(ax.Position(2)-ax.YLim(1)*Dy+Y_points(2)*Dy)-y_begin;
% % Textarrow3.Position=[x_begin y_begin  x_End y_End];