 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 29/05/2020
 function FMMSE=MMSE_MU_MC_MISO(H,P)
        [~,Nt,K,F]=size(H);
        for f=1:F
            Hf=zeros(K,Nt);
            for k=1:K
                Hf(k,:)=H(1,:,k,f);
            end            
%              A=Hf'/(Hf*Hf');
%              alpha=1/trace(A*A');
             Fmmse=  Hf'/(Hf*Hf'+K/P*eye(K));
             Fmmse=sqrt(P)*Fmmse/norm(Fmmse,'fro');
             for k=1:K
                 FMMSE(:,1,k,f)=Fmmse(:,k);
             end
        end
 end