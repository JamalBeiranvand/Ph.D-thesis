clear all
clc
SNR_dB=[-30:2:0];  SNR_linear=10.^(SNR_dB/10.);
N_iter=300; 
Nr=16;Nt=64;
Nt_RF=8;
L=3; % number of rays
M=Nt/Nt_RF; % number of antennas connected to one RF chains
fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;
for i_snr=1:length(SNR_linear)
    i_snr
    SNR=SNR_linear(i_snr);
    temp1=0;temp2=0;temp3=0;temp4=0;temp5=0;
    for i=1:N_iter
        [H,power_matrix,A_BS]=mmWave_channel(Nr,Nt,L,lamada);
%         H=randn(Nr,Nt)+1i*randn(Nr,Nt);
        %%%%%% conventional method %%%%%%%%%%
        F1=SVD_precoding(Nt_RF,H);
        temp1=temp1+log2(det(eye(Nr)+(SNR/8)*H*F1*F1'*H'));
        %%%%% proposed method %%%%%%%%%%%%%%%
        F2=propose_precoding(Nt_RF,Nt,Nr,M,H,SNR);
        temp2=temp2+log2(det(eye(Nr)+(SNR/8)*H*F2*F2'*H'));
        %%%%%% analog method %%%%%%%%%%%%%%%
        F3=full_analog(Nt_RF,Nt,Nr,M,H,SNR);
        temp3=temp3+log2(det(eye(Nr)+(SNR/8)*H*F3*F3'*H'));

        %%%%%% spatially_sparse_precoding %%%%%%%%%%%%%%%%%%%%%%%%%
        [F_RF,F_BB]=spatially_sparse_precoding(Nt_RF,H,A_BS);
        F4=F_RF*F_BB;
        temp4=temp4+log2(det(eye(Nr)+(SNR/8)*H*F4*F4'*H'));
        %%%% SVD %%%%%%%%%%%%%%%%%%%%
        F5=hybrid_precoding(Nt_RF,Nt,Nr,M,H,SNR);
        temp5=temp5+log2(det(eye(Nr)+(SNR/8)*H*F5*F5'*H'));
    end
    C1(i_snr)= real(temp1/N_iter);
    C2(i_snr)= real(temp2/N_iter);     
    C3(i_snr)= real(temp3/N_iter);
    C4(i_snr)= real(temp4/N_iter);
    C5(i_snr)= real(temp5/N_iter);
end
plot(SNR_dB,C1,'r-o','Linewidth',1.5);
hold on
plot(SNR_dB,C2,'b-s','Linewidth',1.5);
% hold on
plot(SNR_dB,C3,'g-^','Linewidth',1.5);
hold on
plot(SNR_dB,C4,'y-^','Linewidth',1.5);
hold on
plot(SNR_dB,C5,'y-^','Linewidth',1.5);
xlabel('SNR (dB)')
ylabel('Capacity')
grid on 