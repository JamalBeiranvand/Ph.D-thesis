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
addpath(('F:\Ph.D\Matlab\Functions') )% Add path of the "Functions" folder.
K=3;                           % number of users
b =2;                    % Values of Bits, it can be a vector, you can use 'inf' for infinity resolution
Txz = 6;                       % Number of the transmitter antennas on z axis
Txy =6;                      % Number of the transmitter antennas on y axis
Rxz = 1;                       % Number of the receiver antennas on z axis
Rxy = 1;                       % Number of the receiver antennas on y axis
% % Ns  = 1;                       % Number of data streams
Nrf = K;                     % Number of RF chains at BS
Ncl =5;                       % Number of Channel clusters(Scatters)
% Nray= 1;                       % Number of rays in each cluster
SNR_dB =0:5:20;             % Signal to noise ratio in dB
TxSector=[90 180];
RxSector=[60 180];
realization = 100;              % Iteration of simulation
%% Initialize Parameters
Tx  = [Txz,Txy];               % Srtucture of the transmit antanna array
Rx  = [Rxz,Rxy];               % Srtucture of the receiver antanna array 
Nt  = Txz*Txy;                 % Number of the transmit antennas
Nr  = Rxz*Rxy;                 % Number of the receive antennas
SNR = 10.^(SNR_dB./10);        
%% Allocate memory space for matrices
R_ZF=zeros(length(SNR_dB),realization);
R_MMSE=zeros(length(SNR_dB),realization);
R_MRC=zeros(K,1);
for Sn=1:length(SNR_dB)
    for it=1:realization
            Phs=1/(sqrt(Nt))*exp(1j*(0:2^b-1)*2*pi/(2^b));
%              [H,At,Ar,~]=MU_Channel(Tx,K,Ncl,'TxSector',TxSector,'RxSector',RxSector);
             [H,~,~,~]=MU_Channel(Tx,K,Ncl); 
             H=H.';
             %% Initialize FRF, FBB.
             q=ones(K,1)/SNR(Sn);
             [~,index]=max(abs(H).'); 
             FRF=zeros(Nt,Nrf);
             for i=1:Nt
             FRF(i,index(i))=1/(sqrt(Nt));
             end
             FBB=sqrt(SNR(Sn))*randn(Nrf,K);

             for TConv=1:20
                 %% Analog beamformer design:
                 SNRMatrix=H'*FRF*FBB;
                 SNRvec=abs(diag(SNRMatrix)).^2;
                 I=sum(abs(SNRMatrix).^2,2)-SNRvec;
                 SINRk=SNRvec./(I+1);
                  R=sum(log2(1+SINRk));
                 for i=1:Nt  
                     R=0;
                     for nrf=1:Nrf
                         for ps=1:length(Phs)
                             FRF(i,:)=0;
                             FRF(i,nrf)=Phs(ps);
                             SNRMatrix=H'*FRF*FBB;
                             SNRvec=abs(diag(SNRMatrix)).^2;
                             I=sum(abs(SNRMatrix).^2,2)-SNRvec;
                             SINRk=SNRvec./(I+1);
                             if sum(log2(1+SINRk))>R
                                R=sum(log2(1+SINRk));
                                nrfs=nrf;
                                Phas=Phs(ps);
                             end             
                         end                        
                     end         
                     FRF(i,:)=0;
                     FRF(i,nrfs)=Phas;                     
                 end
%                  R
                %  sum(FRF~=0,2)

                 %% Digital beamformer design:

                 for k=1:K
                     h_hat(k,:)=H(:,k)'*(FRF);
                 end
                 %% Update uplink beamformer fk as (46) and (47).
                 for Covr=1:30
                     for k=1:K
                         Sig=zeros(K,K);
                         for j=1:K
                             if j~=k
                                 Sig=Sig+q(j)*h_hat(j,:)'*h_hat(j,:);
                             end
                         end
                         fk(:,k)= inv(Sig+eye(K,K))*h_hat(k,:)';
                         fk(:,k)=fk(:,k)/norm(fk(:,k));
                         Eps(k)=h_hat(k,:)*inv(Sig+eye(K,K))*h_hat(k,:)';
%                           Eps(k)=h_hat(k,:)*fk(:,k);
                     end

                     % Update uplink power qk as (48).
                     q=WaterFilling(Eps,SNR(Sn));
                %      sum(q)
                 end
                 %% SINRkUp
                 for k=1:K
                     num=q(k)*fk(:,k)'*h_hat(k,:)'*h_hat(k,:)*fk(:,k);
                     d=zeros(size(h_hat(k,:)'*h_hat(k,:)));
                     for j=1:K
                         if j~=k
                             d=d+q(j)*h_hat(j,:)'*h_hat(j,:);
                         end
                     end
                     dem=fk(:,k)'*d*fk(:,k)+1;
                     SINRkUp(k)=num/dem;
                 end
                 SINRz=find(SINRkUp==0);
                 SINRnz=find(SINRkUp~=0);
                %% Compute downlink beamformer FBB by (49)-(51).
                Bmatix=zeros(K,K);
                for i=1:K
                    for j=1:K
                        if i==j
                            if SINRkUp(i)==0
%                                 x=1
                                B(i,j)=0;
                            else
                                B(i,j)=abs(h_hat(i,:)*fk(:,j))/SINRkUp(i);
                            end
                        else
                                B(i,j)=-abs(h_hat(i,:)*fk(:,j))^2;
                        end
                    end
                end
                 
%                 for z=1:length(SINRz)
%                     B(z,:)=[];
%                     B(:,z)=[];
%                 end
                Binv=inv(B.');                
                for k=K
                    pk(k)=sum(Binv(k,:));
                end
                pk(SINRz)=0;
                find(pk==0)
                for k=1:K
                    FBBhat(:,k)=sqrt(pk(k))*fk(:,k);
                end
                %% Normalize FBB as (52).
%                 FBB=sqrt(SNR(Sn)*K)*FBBhat/norm(FRF*FBBhat,'fro');
                  FBB=sqrt(SNR(Sn))*FBBhat/norm(FRF*FBBhat,'fro');
             end
%     SNR(Sn)        
%     sum(sum(abs(FRF*FBB).^2))
    SNRMatrix=H'*FRF*FBB;
    SNRvec=abs(diag(SNRMatrix)).^2;
    I=sum(abs(SNRMatrix).^2,2)-SNRvec;
    SINRk=SNRvec./(I+1);
    RR(it)=sum(log2(1+SINRk));
     end
     Rsum(Sn)=mean(RR)
end
plot(SNR_dB,Rsum)
axis([min(SNR_dB) max(SNR_dB) 5 30])