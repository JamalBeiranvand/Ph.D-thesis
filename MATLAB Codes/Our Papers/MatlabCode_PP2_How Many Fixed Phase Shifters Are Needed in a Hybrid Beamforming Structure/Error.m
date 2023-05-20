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
K=20;                           % number of users
NPS=17;
Txz = 12;                       % Number of the transmitter antennas on z axis
Txy =12;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 1;                       % Number of the receiver antennas on y axis
% Ns  = 1;                       % Number of data streams
% NrfTx = 9;                     % Number of RF chains at BS
Ncl =10;                       % Number of Channel clusters(Scatters)
Nray= 5;                       % Number of rays in each cluster
AngSpread = 10;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
TxSector=[360 360];
RxSector=[360 360];
SNR_dB = 0;             % Signal to noise ratio in dB
realization = 1000;              % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10); 
%%
[Setb1,Sb1]=GenUniqueSet(NPS(1));
% [Setb2,Sb2]=GenUniqueSet(NPS(2));
% [Setb3,Sb3]=GenUniqueSet(NPS(3));
% [Setb4,Sb4]=GenUniqueSet(NPS(4));
for reali=1:realization
        % generate channel
        [H,At,Ar,~]=MU_Channel(Tx,K,Ncl,'Nray',Nray,'AngSpread',AngSpread);
        G= diag(ones(K,1))*sqrt(SNR);   
        % Full Digital ZF beamforming
        Z=ZF_Down(H);
%         *max(abs(Setb1))/max(max(abs(Z)))
         Z=Z*max(abs(Setb1))/max(max(abs(Z)));
        Fzf(:,reali)=reshape(Z,[],1);        
         % Implement 
        [FRF_ZF,~]=MUMISO_DownModified(Z,Setb1,Sb1);
%         FRF_ZF=FRF_ZF/sqrt(sum(sum(abs(FRF_ZF).^2))
        E_ZF(:,reali)=reshape(FRF_ZF,[],1);
        
        E(:,reali)=reshape(abs(Z-FRF_ZF),[],1);
%         hist(abs(reshape(Fzf,1,[])),100)
%         hold on
end
% Fzf=reshape(Fzf,1,[]);
% hist(abs(Fzf),100)
% hold on
% figure
% E_ZF=reshape(E_ZF,1,[]);
% hist(abs(E_ZF),40)
%  figure
%  hist(abs(Setb1),10000)
 
 figure
 E=reshape(E,1,[]);
%  mean(E)
 mean(E)/max(abs(Setb1))
 hist(abs(E),100)
 
