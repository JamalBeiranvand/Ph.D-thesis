clear
clc
addpath(genpath('D:\Data Set\Data1' ))
filename = 'DNN1.h5';
lgraph = importKerasLayers(filename,'ImportWeights',true);
net = assembleNetwork(lgraph);
load("TestingData.mat")
load("TestingSwitchMatrix.mat")
realization = size(TestingData,1);              % Iteration of simulation
Nt  = size(TestingData,2);                 % Number of the transmit antennas
K=size(TestingData,3); % number of users
Ncom=sum(sum(reshape(TestingSwitchMatrix(1,:,:),Nt,[]),2)==K); 
NAntenna_User=sum(reshape(TestingSwitchMatrix(1,:,:),Nt,[]))-Ncom;

Error1=zeros(realization,K);

for i=1:realization
    TrueLabes=reshape(TestingSwitchMatrix(i,:),Nt,[]);
    NN_Probability=predict(net,squeeze(TestingData(i,:,:,:)));
    NN_Probability=reshape(NN_Probability,Nt,[]);
    [~,~,NN_pre]=Predict_NN(NN_Probability,Ncom,NAntenna_User);
    Error1(i,:)=sum(NN_pre~=TrueLabes);
end

% 


% 
% Ncom=sum(sum(squeeze(TestingSwitchMatrix(1,:,:)),2)==K); 
% NAntenna_User=sum(reshape(TestingUserLabels(1,:),Nt,[]));