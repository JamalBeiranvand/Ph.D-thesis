 
% Titel: "Near-Optimal Hybrid Analog and Digital Precoding for Downlink mmWave Massive MIMO Systems"
% Address: https://ieeexplore.ieee.org/document/7248508
% Input: 
% Number of data streams = Number of RF chains Ns;
% SNR  
% Channel matrix H
function [FRF,FBB ] = PartSIC(H,Ns,SNR,Nq)
[~,Nt] = size(H);
M = Nt/Ns;
A = [];
D = [];
PSIC = zeros(M*Ns,Ns);
% G = R*(H)'*H*(R)';
for n = 1:Ns
    G=[];
    Tm = eye(size(H,1))+ SNR/Ns * H*PSIC(:,1:(n-1))*PSIC(:,1:(n-1))'*H';
    G=H'*inv(Tm)*H;
    SubChannel=G(((n-1)*M)+1:n*M,((n-1)*M)+1:n*M);
    [~,~,V] = svd(SubChannel);
    v1 = V(:,1);    
    A = blkdiag( A, exp(1i*angle(v1)) );
    D = blkdiag( D, real(v1'*exp(1i*angle(v1)))/(M) );
    PSIC(:,n)=[zeros(M*(n-1),1);D(n,n)*exp(1i*angle(v1));zeros(M*(Ns-n),1)];
    Popt(:,n)=[zeros(M*(n-1),1); v1; zeros(M*(Ns-n),1)];
end
FRF=A;
FBB=D;
  %% quantization
    if (exist('Nq','var'))
        FRF=Qtz(FRF,Nq);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    end


