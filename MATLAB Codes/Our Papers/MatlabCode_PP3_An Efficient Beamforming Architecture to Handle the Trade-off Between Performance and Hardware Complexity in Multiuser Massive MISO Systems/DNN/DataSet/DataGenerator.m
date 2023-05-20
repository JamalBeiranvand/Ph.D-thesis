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
 addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Our Papers\Matlabcodes__Dynamic Structure' ))% Add path of the "Functions" folder.
 addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' ))% Add path of the "Functions" folder.

K=10; % number of users
Ncom=113;
LSNR=[15,20,25];
LIteration = 1000;              % Iteration of simulation
LC=[5];
LN=33; % noisy version for each level of noise


Txz = 16;                       % Number of the transmitter antennas on z axis
Txy =16;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 1;                       % Number of the receiver antennas on y axis
% Ns  = 1;                       % Number of data streams
% NrfTx = 9;                     % Number of RF chains at BS
%Ncl =5;                       % Number of Channel clusters(Scatters)
Nray= 10;                       % Number of rays in each cluster
AngSpread = 10;       % Angle spread, standard deviation of the angles in azimuth and elevation both of Rx and Tx
TxSector=[60 360];
RxSector=[20 360];
% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
MultiLabels=zeros([(LN*length(LSNR)+1)*length(LC)*LIteration,Nt]);
UserLabels=zeros([(LN*length(LSNR)+1)*length(LC)*LIteration,Nt,K]);
SwitchMatrix=zeros([(LN*length(LSNR)+1)*length(LC)*LIteration,Nt,K]);
Ht=zeros([(LN*length(LSNR)+1)*length(LC)*LIteration,Nt,K,3]);

indx=0;
for N=1:LIteration
    for lc=1:length(LC)
        Ncl=LC(lc);
        [H,At,Ar,~]=MU_Channel(Tx,K,Ncl,'Nray',Nray,'AngSpread',AngSpread); 
        H=H.';
        [ComSet,~]=ComAntSelection_Greedy(H,Ncom);        
        KSets=nComAntSelection_Greedy(H,ComSet);

        for ln=1:LN+1
             if ln==1 
                 indx=indx+1;
                 Ht(indx,:,:,1)=abs(H);
                 Ht(indx,:,:,2)=real(H);
                 Ht(indx,:,:,3)=imag(H);
                MultiLabels(indx,ComSet)=1;
                for k=1:K
                    UserLabels(indx,KSets(KSets(:,k)>0,k),k)=1;
                    SwitchMatrix(indx,KSets(KSets(:,k)>0,k),k)=1;
                end
                SwitchMatrix(indx,ComSet,:)=1;
             else
                 for s=1:length(LSNR)
                     indx=indx+1;
                     snr=LSNR(s);
                     HN=awgn(H,snr);
                     Ht(indx,:,:,1)=abs(HN);
                     Ht(indx,:,:,2)=real(HN);
                     Ht(indx,:,:,3)=imag(HN); 

                    MultiLabels(indx,ComSet)=1;
                    for k=1:K
                        UserLabels(indx,KSets(KSets(:,k)>0,k),k)=1;
                        SwitchMatrix(indx,KSets(KSets(:,k)>0,k),k)=1;
                    end
                    SwitchMatrix(indx,ComSet,:)=1;
                 end
             end
        end
    end
end
Ht=single(Ht);
MultiLabels=single(MultiLabels);
UserLabels=single(UserLabels);
SwitchMatrix=single(SwitchMatrix);

TestIndex=1:indx;
TraningIndex=randperm(indx,floor(2*indx/3));
TestIndex(TraningIndex)=[];

TrainingData=Ht(TraningIndex,:,:,:);
TrainingSwitchMatrix=SwitchMatrix(TraningIndex,:,:);
TrainingSwitchMatrix=reshape(TrainingSwitchMatrix,[length(TraningIndex),Nt*K]);

TrainingLabel=MultiLabels(TraningIndex,:);
TrainingUserLabels=UserLabels(TraningIndex,:,:);
TrainingUserLabels=reshape(TrainingUserLabels,[length(TraningIndex),Nt*K]);


TestingData=Ht(TestIndex,:,:,:);
TestingSwitchMatrix=SwitchMatrix(TestIndex,:,:);
TestingSwitchMatrix=reshape(TestingSwitchMatrix,[length(TestIndex),Nt*K]);
TestingLabel=MultiLabels(TestIndex,:);
TestingUserLabels=UserLabels(TestIndex,:,:);
TestingUserLabels=reshape(TestingUserLabels,[length(TestIndex),Nt*K]);


save("TrainingData.mat","TrainingData")
save("TrainingSwitchMatrix.mat","TrainingSwitchMatrix")
save("TrainingLabel","TrainingLabel")
save("TrainingUserLabels","TrainingUserLabels")

save("TestingData","TestingData")
save("TestingSwitchMatrix.mat","TestingSwitchMatrix")
save("TestingLabel","TestingLabel")
save("TestingUserLabels","TestingUserLabels")



