function [ WBB,WRF,FRF,FBB] = PE_AltMin_OFDM_MIMO( H, NRF,Ns,P )

[Nr, Nt, F] = size(H);
    for f = 1:F
        [Uj,~,Vj] = svd(H(:,:,f));
        Fopt(:,:,f) = Vj(1:Nt,1:Ns);
        Wopt(:,:,f) = Uj(1:Nr,1:Ns);
    end


mynorm = [];
FRF = exp( 1i * unifrnd (0,2*pi,Nt,NRF) );
% FBB = zeros(NRF, Ns, K);
while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-1)
    mynorm = [0,0];
    temp = zeros(Nt, NRF);
    for f = 1:F
        [U,~,V] = svd(Fopt(:,:,f)'*FRF);
        FBB(:,:,f) = V(:,[1:Ns])*U';
        mynorm(1) = mynorm(1) + norm(Fopt(:,:,f) * FBB(:,:,f)' - FRF,'fro')^2;
        temp = temp + Fopt(:,:,f) * FBB(:,:,f)';
    end

    FRF = exp(1i * angle(temp));
    for f = 1:F
        mynorm(2) = mynorm(2) + norm(Fopt(:,:,f) * FBB(:,:,f)' - FRF,'fro')^2;
    end
end


    for f = 1:F
        FBB(:,:,f) = sqrt(P) * FBB(:,:,f) / norm(FRF * FBB(:,:,f),'fro');
    end
    
    %%
    
mynorm = [];
WRF = exp( 1i * unifrnd (0,2*pi,Nr,NRF) );
% FBB = zeros(NRF, Ns, K);
while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-1)
    mynorm = [0,0];
    temp = zeros(Nr, NRF);
    for f = 1:F
        [U,~,V] = svd(Wopt(:,:,f)'*WRF);
        WBB(:,:,f) = V(:,[1:Ns])*U';
        mynorm(1) = mynorm(1) + norm(Wopt(:,:,f) * WBB(:,:,f)' - WRF,'fro')^2;
        temp = temp + Wopt(:,:,f) * WBB(:,:,f)';
    end

    WRF = exp(1i * angle(temp));
    for f = 1:F
        mynorm(2) = mynorm(2) + norm(Wopt(:,:,f) * WBB(:,:,f)' - WRF,'fro')^2;
    end
end
    
end