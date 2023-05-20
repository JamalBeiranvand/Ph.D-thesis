% This function is written according to Algorithm2 on page 7 in:
% "Spatially Sparse Precoding in Millimeter Wave MIMO Systems" 
% you can see the paper on this link: 
% https://ieeexplore.ieee.org/document/6717211
% Inputs:
    % Fopt:  the optimal fully digital precoder matrix
    % NRF: Nuber of RF chains
    % At: transmit array response vectors
% outputs
    % WRF: RF combiner matrix
    % WBB: baseband combiner matrix
%%
function [WRF,WBB,FRF,FBB] = OMP_Alg(H,At,Ar,NRF,Ns,SNR,Nq)
    NtRF=NRF;
    NrRF=NRF;
    noisevar=1/SNR;
    [~,~,v] = svd(H);
    Fopt = v(:,1:Ns);
    [Nr,Nt] = size(H);

    Frf = complex(zeros(Nt,NtRF));
    Fres = Fopt;
    for m = 1:NtRF
        Psi = At'*Fres;
        [~,k] = max(diag(Psi*Psi'));
        Frf(:,m) = At(:,k);
        Fbb = (Frf(:,1:m)'*Frf(:,1:m))\Frf(:,1:m)'*Fopt;
        temp = Fopt-Frf(:,1:m)*Fbb;
        Fres = temp/norm(temp,'fro');
    end
    Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');

    Wmmse = ((Fbb'*Frf'*(H'*H)*Frf*Fbb+noisevar*Ns*eye(Ns))\Fbb'*Frf'*H')';
    %Wmmse = (1/Ns*(Fbb'*Frf'*H')/(1/Ns*H*(Frf*(Fbb*Fbb')*Frf')*H'+noisevar*eye(Nr)))';
    Wrf = complex(zeros(Nr,NrRF));
    Wres = Wmmse;
    Ess = 1/Ns*eye(Ns);
    Eyy = H*Frf*Fbb*Ess*Fbb'*Frf'*H'+noisevar*eye(Nr);
    for m = 1:NrRF
        Psi = Ar'*Eyy*Wres;
        [~,k] = max(diag(Psi*Psi'));
        Wrf(:,m) = Ar(:,k);
        Wbb = (Wrf(:,1:m)'*Eyy*Wrf(:,1:m))\(Wrf(:,1:m)'*Eyy*Wmmse);
        temp = Wmmse-Wrf(:,1:m)*Wbb;
        Wres = temp/norm(temp,'fro');
    end

    FRF=Frf;
    FBB=Fbb;
    WRF=Wrf;
    WBB=Wbb;
    
        % quantization
    if (exist('Nq','var'))
        FRF=Qtz(FRF,Nq);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        WRF=Qtz(WRF,Nq);
    end
end

