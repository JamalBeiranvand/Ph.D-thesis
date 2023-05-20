% This function is written according to PE-AltMin Algorithm on page 8 in:
    % Titel: "Alternating Minimization Algorithms for Hybrid Precoding in Millimeter Wave MIMO Systems"
    % Address: https://ieeexplore.ieee.org/document/7397861
% Inputs:
    % Fopt:  the optimal fully digital precoder matrix
    % NRF: Nuber of RF chains
% outputs
    % FRF: RF combiner matrix
    % FBB: baseband combiner matrix
function [WRF,WBB,FRF,FBB]= PE_AltMin(H,NRF,Ns,Nq)
        [U,~,V] = svd(H);
        Fopt   = V(:,1:Ns);  
        Wopt   = U(:,1:Ns);
        [Nr,Nt] = size(H);
        
    mynorm = [];
    %Construct FRF with random phases and set k = 0;
    FRF = exp( 1i * unifrnd (0,2*pi,Nt,NRF) );
    % stopping criterion: abs(mynorm(end)-mynorm(end-1) > 1e-3
    while (isempty(mynorm) || abs( mynorm(end) - mynorm(end-1) ) > 1e-3)
        % compute the SVD        
        [U,~,V] = svd(Fopt'*FRF);
        % Compute FBB according to eq (28) in page 8
        FBB = V(:,1:Ns)*U';
        % stopping criterion
        mynorm = [mynorm, norm(Fopt * FBB' - FRF,'fro')^2];
        % extracted the phases of FRF  from the phases of an equivalent precoder
        FRF = exp(1i * angle(Fopt * FBB'));
        % stopping criterion
        mynorm = [mynorm, norm(Fopt * FBB' - FRF,'fro')^2];
    end
    FBB = sqrt(Ns)*FBB/norm(FRF*FBB,'fro');
  %% receiver  
    mynorm = [];
    U=[];
    V=[];
    %Construct FRF with random phases and set k = 0;
    WRF = exp( 1i * unifrnd (0,2*pi,Nr,NRF) );
    % stopping criterion: abs(mynorm(end)-mynorm(end-1) > 1e-3
    while (isempty(mynorm) || abs( mynorm(end) - mynorm(end-1) ) > 1e-3)
        % compute the SVD        
        [U,~,V] = svd(Wopt'*WRF);
        % Compute FBB according to eq (28) in page 8
        WBB = V(:,1:Ns)*U';
        % stopping criterion
        mynorm = [mynorm, norm(Wopt * WBB' - WRF,'fro')^2];
        % extracted the phases of FRF  from the phases of an equivalent precoder
        WRF = exp(1i * angle(Wopt * WBB'));
        % stopping criterion
        mynorm = [mynorm, norm(Wopt * WBB' - WRF,'fro')^2];
    end
    %% quantization
    if (exist('Nq','var'))
        FRF=Qtz(FRF,Nq);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        WRF=Qtz(WRF,Nq);
    end
    
end