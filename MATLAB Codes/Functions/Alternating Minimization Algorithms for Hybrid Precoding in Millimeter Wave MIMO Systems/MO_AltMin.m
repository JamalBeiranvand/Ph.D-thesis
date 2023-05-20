
% MO-AltMin Algorithm: Manifold Optimization Based Hybrid Precoding for the Fully-connected Structure
% this functions is written according "MO-AltMin Algorithm" on page 6 of 
    % Titel: "Alternating Minimization Algorithms for Hybrid Precoding in Millimeter Wave MIMO Systems"
    % Address: https://ieeexplore.ieee.org/document/7397861
 % you need to "sig_manif" function and manopt package   
function [WRF,WBB, FRF,FBB ] = MO_AltMin( H, NRF,Ns,Nq )
        [U,~,V] = svd(H);
        Fopt = V(:,1:Ns);  
        Wopt = U(:,1:Ns);
    %% FRF,FBB    
    [Nr, Nt] = size(H);
    y = [];
    FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    while(isempty(y) || abs(y(1)-y(2))>1e-3)
            FBB = pinv(FRF) * Fopt;
            y(1) = norm(Fopt - FRF * FBB,'fro')^2;
            [FRF, y(2)] = sig_manif(Fopt, FRF, FBB);
    end
    
    %% WRF,WBB
    y = [];
    WRF = exp( 1i*unifrnd(0,2*pi,Nr,NRF) );
    while(isempty(y) || abs(y(1)-y(2))>1e-3)
            WBB = pinv(WRF) * Wopt;
            y(1) = norm(Wopt - WRF * WBB,'fro')^2;
            [WRF, y(2)] = sig_manif(Wopt, WRF, WBB);
    end
      %% quantization
    if (exist('Nq','var'))
        FRF=Qtz(FRF,Nq);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        WRF=Qtz(WRF,Nq);
    end
end