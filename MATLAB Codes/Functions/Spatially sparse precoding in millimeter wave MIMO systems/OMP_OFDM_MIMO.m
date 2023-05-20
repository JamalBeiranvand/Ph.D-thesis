function [ WBB,WRF, FRF,FBB ] = OMP_OFDM_MIMO( H, NRF,Ns,P, At,Ar )

[Nr, Nt, F] = size(H);
    for f = 1:F
        [Uj,~,Vj] = svd(H(:,:,f));
        Fopt(:,:,f) = Vj(1:Nt,1:Ns);
        Wopt(:,:,f) = Uj(1:Nr,1:Ns);
    end
    
FRF = [];
% FBB=[];
Fres = Fopt;
for i = 1:NRF
    temp = 0;
    for f = 1:F
        PU(:,:,f) = At' * Fres(:,:,f);
        temp = temp + sum( abs(PU(:,:,f)).^2, 2 );
    end
    [aa,bb] = max(temp);
    FRF = [FRF , At(:,bb)];
    for f = 1:F
        FB{f} = pinv(FRF) * Fopt(:,:,f);
        Fres(:,:,f) = (Fopt(:,:,f) - FRF * FB{f}) / norm(Fopt(:,:,f) - FRF * FB{f},'fro');
    end
end

 for f = 1:F
      FBB(:,:,f) = sqrt(P) * FB{f} / norm(FRF * FB{f},'fro');
 end
 
 %%
 WRF = [];
Fres = Wopt;
for i = 1:NRF
    temp = 0;
    for f = 1:F
        PU(:,:,f) = Ar' * Fres(:,:,f);
        temp = temp + sum( abs(PU(:,:,f)).^2, 2 );
    end
    [aa,bb] = max(temp);
    WRF = [WRF , Ar(:,bb)];
    for f = 1:F
        WB {f}= pinv(WRF) * Wopt(:,:,f);
        Fres(:,:,f) = (Wopt(:,:,f) - WRF * WB{f}) / norm(Wopt(:,:,f) - WRF * WB{f},'fro');
    end
end
 for f = 1:F
      WBB(:,:,f) = sqrt(P) * WB{f} / norm(FRF * WB{f},'fro');
 end

 

end