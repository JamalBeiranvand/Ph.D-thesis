function R=SumRate_MIMO_OFDM_Hybrid(H,WBB,WRF,FRF,FBB,G)
[~,~,F]=size(H);
[~,Ns,~]=size(WBB);
R=0;
      for f=1:F  
          Wf=WRF*WBB(:,:,f);
          Ff=FRF*FBB(:,:,f);
          Hf= H(:,:,f);             
          R= R+log2(det(eye(Ns) + G* pinv(Wf)  * Hf * Ff  * (Ff)'  * Hf' * Wf))/(F);        
      end
end