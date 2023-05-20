function FZF=ZF_MU_MC_MISO(H,P)
        [~,Nt,K,F]=size(H);
        for f=1:F
            Hf=zeros(K,Nt);
            for k=1:K
                Hf(k,:)=H(1,:,k,f);
            end            
%              A=Hf'/(Hf*Hf');
%              alpha=1/trace(A*A');
             Fzf=  Hf'/(Hf*Hf');
             Fzf=sqrt(P)*Fzf/norm(Fzf,'fro');
             for k=1:K
                 FZF(:,1,k,f)=Fzf(:,k);
             end
        end
end