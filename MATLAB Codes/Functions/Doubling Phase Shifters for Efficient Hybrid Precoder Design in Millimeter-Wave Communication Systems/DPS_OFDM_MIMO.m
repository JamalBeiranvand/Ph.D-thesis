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
function [WBB,WRF,FRF,FBB]=DPS_OFDM_MIMO(H,NRF,Ns,P)
    [Mr,Mt,F]=size(H);
    FOP=[];
    WOP=[];
    for f=1:F
        [Wop,~,Fop]=svd(H(:,:,f));
        FOP=[FOP,Fop(:,1:Ns)];
        WOP=[WOP,Wop(:,1:Ns)];
    end
    
    [UF,S,VF]=svd(FOP); 
    FRF=UF(:,1:NRF);
    Fbb=S(1:NRF,1:NRF)*VF(:,1:NRF)';
    for f=1:F
        FBB(:,:,f)=Fbb(:,(f-1)*Ns+1:f*Ns);  
        FBB(:,:,f)=sqrt(P)*FBB(:,:,f)/norm(FRF*FBB(:,:,f),'fro');
    end
    
    [UF,S,VF]=svd(WOP);
    WRF=UF(:,1:NRF);
    Wbb=S(1:NRF,1:NRF)*VF(:,1:NRF)';
   for f=1:F
        WBB(:,:,f)=Wbb(:,(f-1)*Ns+1:f*Ns);  
   end

end