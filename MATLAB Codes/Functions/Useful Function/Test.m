% this function plot the Figure 3 on page 12 of     
    % Titel: "Alternating Minimization Algorithms for Hybrid Precoding in Millimeter Wave MIMO Systems"
    % Address: https://ieeexplore.ieee.org/document/7397861
% you need below functions and packages to run this code
    % ChannelUSPA & laprnd: to generate channel
    % OMP_Precoder & OMP_Combiner: 
    % Anolog_Beamforming: 
    % SDR_AltMin & 
    % SIC
    % MO_AltMin , sig_manif & manopt (package)
    
%  NOTE: For the "SDR_AltMin" and "SIC" algorithms, the authors used a cvx (Version 2.1, Build 1103 (9714d49)). Users should include a cvx package when these algorithms are implemented.   
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
Bits = 3;                      % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 8;                      % Number of the transmitter antennas on z axis
Txy = 8;                      % Number of the transmitter antennas on y axis
Rxz = 4;                       % Number of the receiver antennas on z axis
Rxy = 4;                       % Number of the receiver antennas on y axis
Ns  = 1;                       % Number of data streams
NRF = 4;                      % Number of RF chains
Ncl =8;                       % Number of Channel clusters(Scatters)
Nray= 10;                       % Number of rays in each cluster
AngSpread = 7.5/180*pi;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
SNR_dB = -40:5:0;             % Signal to noise ratio in dB
realization = 500;             % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
Fopt=zeros(Nt,Ns);                          % Fopt:  the optimal fully digital precoder matrix,  
Wopt=zeros(Nr,Ns);                          % Wopt: the optimal fully digital combiner
R_OMP=zeros(length(SNR_dB),realization);    % allocate space to save Spectral efficiency of OMP algorithm
R_opt=zeros(length(SNR_dB),realization);    % allocate space to save Spectral efficiency of Optimum 
R_Ang=zeros(length(SNR_dB),realization);    % allocate space to save Spectral efficiency of Analog Beamforming
R_SIC=zeros(length(SNR_dB),realization);    % allocate space to save Spectral efficiency of SIC Method
R_SDR=zeros(length(SNR_dB),realization);    % allocate space to save Spectral efficiency of SDR Method
R_MO=zeros(length(SNR_dB),realization);     % allocate space to save Spectral efficiency of SDR Method
%% Main code 
h = waitbar(0,'Please wait...');% for display Time bar
tic
for r = 1:length(SNR)
    for reali = 1:realization
        % generate channel
         [H,At,Ar,alpha]=Channel(Tx,Rx,Ncl,'Nray',Nray);
         
        %Full Digital beamforming
        [Wopt,Fopt]=DigitalBeamforming(H,Ns);
        
        % obtain FRF and FBB by function of OMP algorithm
        [WRF_OPM,WBB_OPM,FRF_OMP, FBB_OMP] = OMP_Alg(H,At,Ar,NRF,Ns,SNR(r));
        
%         % obtain FRF and WRF by function of Analog Beamforming
%         [WRF_Ang,FRF_Ang] = AnalogBeamforming(At,Ar,alpha,Ns);
%         
%         % obtain FRF and FBB by function of SDR_AltMin
%        [FRF_SDR, FBB_SDR] = PartSDR_AltMin(H,NRF,Ns);
%        
%         % obtain FRF and WRF by function of SIC
%         [FRF_SIC, FBB_SIC] = PartSIC( H, Ns, SNR(r));
%         
% %         % obtain FRF and FBB by function of MO_AltMi 
%             [WRF_MO, WBB_MO,FRF_MO, FBB_MO]=MO_AltMin(H,NRF,Ns);
%         [FRF_MO, FBB_MO] = MO_AltMin( Fopt, NRF);
%         FBB_MO = sqrt(Ns) * FBB_MO / norm(FRF_MO * FBB_MO,'fro');
%         % obtain WRF and WBB by function of MO_AltMi 
%         [WRF_MO, WBB_MO] = MO_AltMin(Wopt, NRF);
        
         % compute Spectral efficiency for each algorithm 
         R_OMP(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_OPM * WBB_OPM)* H * FRF_OMP * (FBB_OMP) * FBB_OMP' * FRF_OMP' * H' * WRF_OPM * WBB_OPM));
%          R_MO(r,reali)  = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_MO  * WBB_MO) * H * FRF_MO * FBB_MO * (FBB_MO)' * FRF_MO' * H' * WRF_MO * WBB_MO));
%          R_Ang(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(WRF_Ang)* H * FRF_Ang * (FRF_Ang)' * H' * WRF_Ang) );
%          R_SDR(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(Wopt)   * H * FRF_SDR * FBB_SDR * (FBB_SDR)' * FRF_SDR' * H' * Wopt));
%          R_SIC(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(Wopt)   * H * FRF_SIC * FBB_SIC * (FBB_SIC)' * FRF_SIC' * H' * Wopt));
        % compute Optimal Spectral efficiency  
         R_opt(r,reali) = log2(det(eye(Ns) + SNR(r)/Ns * pinv(Wopt) * H * Fopt * (Fopt)' * H' * Wopt));
         
         % for display Time bar 
         h = waitbar(((r-1)*realization+reali)/(length(SNR)*realization),h,['%',num2str(((r-1)*realization+reali)/(length(SNR)*realization)*100,'%4.1f')]);

    end
end
toc
%% Plot results 
figure 
grid on
hold on
plot(SNR_dB,sum(R_opt,2)/realization,'r-o','MarkerSize',10,'LineWidth',2)
plot(SNR_dB,sum(R_OMP,2)/realization,'Marker','^','MarkerSize',10,'LineWidth',2,'Color',[0 .5 0])
%  plot(SNR_dB,sum(R_SDR,2)/realization,'Marker','diamond','MarkerSize',10,'LineWidth',2,'Color',[0.8 0.5 0]);
% plot(SNR_dB,sum(R_SIC,2)/realization,'c-+','MarkerSize',10,'LineWidth',2,'Color',[.1 .2 .5]);
% plot(SNR_dB,sum(R_Ang,2)/realization,'b-*','LineWidth',2,'MarkerSize',10,'Color',[.3 .2 .5])
%  plot(SNR_dB,sum(R_MO,2)/realization,'r-p','MarkerSize',10,'LineWidth',2,'Color',[.7 0 1])

xlabel('SNR (dB)','FontSize',20,'FontWeight','normal','Color','k')
ylabel ('Espectral Efficiency (bits/sec/Hz)','FontSize',20,'FontWeight','normal','Color','k')
% legend({'Optimal Digital Precoder','OMP Algorithm','SDR','SIC','Analog Beamforming'},'Location','southeast','FontSize',25)
ax = gca; % current axes
ax.FontSize = 20;
axis([min(SNR_dB) max(SNR_dB) 0 max(real(sum(R_opt,2))/realization)])

