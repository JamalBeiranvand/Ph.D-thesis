function [F,HD,Fn,Fc]=SPD(H,A,SetK)
  [Nt,K]=size(H);
    if Nt<K
        error('size of the channel matrix is not true')
    end
    Nuk=sum(SetK~=0);
    Nu=sum(Nuk);
    HD=(A).'*H;
    %% comput Fn
%     Hn=HD(1:Nu,:);
%     Fn=zeros(sum(Nuk),K);
%     nt=1;
%     for k=1:K    
%         for nuk=1:Nuk(k)
%             Fn(nt,k)=conj(Hn(nt,k));
%             nt=nt+1;
%         end
%     end
%     
%     %% Compute Fc    
%     G=Hn.'*Fn;
%     Hc=HD(Nu+1:end,:);
%     HcDiag=-diag(Hc'*Hc);
%     G=(G-diag(diag(G)));
%     G=G+diag(HcDiag);
%     Fc=-conj(Hc)*((Hc.'*conj(Hc))\G);
%     %trace(Fc'*Fc)
%     PFc=trace(G*G'*(inv(Hc.'*conj(Hc)))');
%     %% Beamforming Matrix F
%     F=[Fn;Fc];
%     F=F/sqrt(sum(sum(abs(F).^2)));
%%
KSets_Fixed=zeros(size(SetK));
KSets_Fixed(SetK~=0)=1:Nu;
[F,Fn,Fc]=SPF(HD,KSets_Fixed);
end
    
    
    