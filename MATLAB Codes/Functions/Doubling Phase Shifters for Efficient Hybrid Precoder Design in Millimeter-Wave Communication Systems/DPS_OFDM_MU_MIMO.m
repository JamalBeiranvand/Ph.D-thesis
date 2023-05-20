% This function is written according to
% "Doubling Phase Shifters for Efficient Hybrid Precoder Design in Millimeter-WaveCommunication Systems" 
% you can see the paper on this link: 
% 
% Inputs:
    % H:  Channel matrix
    % Ns=  Number of streams 
    % NRF: Nuber of RF chains
    % P: Power 
    % b: the number of bits in the resolution of phase shifters
% outputs
    % WRF: RF combiner matrix
    % WBB: baseband combiner matrix
    % FRF: RF precoder matrix
    % FBB: baseband precoder matrix
    
% Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
% google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
%%
function [WBB,WRF,FRF,FBC,FBB]=DPS_OFDM_MU_MIMO(H,NRFTx,NRFRx,Ns,P)
    [Mr,Mt,K,F]=size(H);
    [Wopt,Fopt]=BD(H,Ns,P/K);
    %[Mr,Ns,K,F]=size(Wopt);
    %[Mt,KNs,K,F]=size(Fopt);
    FOP=[];
    WOP=[];
    for k=1:K
        Wk=[];
        for f=1:F
            FOP=[FOP,Fopt(:,:,k,f)];
            Wk=[Wk,Wopt(:,:,k,f)];
        end
        WOP(:,:,k)=Wk;
    end
    
    [UF,S,VF]=svd(FOP); 
    FRF=UF(:,1:NRFTx);
    Fbb=S(1:NRFTx,1:NRFTx)*VF(:,1:NRFTx)';
    for k=1:K
        for f=1:F
            SColumn=(k-1)*F*Ns+(f-1)*Ns+1;
            EColumn=(k-1)*F*Ns+f*Ns;
            FBB(:,:,k,f)=Fbb(:,SColumn:EColumn);  
            FBB(:,:,k,f)=sqrt(Ns*P)*FBB(:,:,k,f)/norm(FRF*FBB(:,:,k,f),'fro');
        end
    end
    
       
   for k=1:K  
       [UF,S,VF]=svd(WOP(:,:,k));        
       Wbb=S(1:NRFRx,1:NRFRx)*VF(:,1:NRFRx)';
       WRF(:,:,k)=UF(:,1:NRFRx);
       for f=1:F
            SColumn=(f-1)*Ns+1;
            EColumn=f*Ns;
            WBB(:,:,k,f)=Wbb(:,SColumn:EColumn);  
       end
   end
%    for f=1:f
%         for k=1:K
%             [W,S,~]=svd(H(:,:,k,f));
%             WRF(:,:,k,f)=W(:,1:Ns);
%             WBB(:,:,k,f)=eye(Ns);
%         end
%    end
   
       FBC= BD_Hybrid(H,WBB,WRF,FRF,P);

end