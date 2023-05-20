function R=SumRate_MUMIMO_Hybrid(H,WBB,WRF,FRF,FBB,G)
[Nr,~,K,F]=size(H);
[~,Ns,~,~]=size(WBB);
R=0;
      for f=1:F  
          for k=1:K
              Wkf=WRF(:,:,k)*WBB(:,:,k,f);
              Fkf=FRF*FBB(:,:,k,f);
              Hkf= H(:,:,k,f);
              
              OMEGA=zeros(Nr);
              for kp=1:K
                  Fkp=FRF*FBB(:,:,kp,f)*G;
                  OMEGA=OMEGA+Hkf*Fkp*(Fkp')*Hkf';
              end               
              Omega=OMEGA-Hkf*Fkf*(Fkf')*Hkf';
              Omega=Wkf'*(Omega+eye(Nr))*Wkf;              
              R= R+log2(det(eye(Ns) + G* Wkf'  * Hkf * Fkf  * (Fkf)'  * Hkf' * Wkf/(Omega)))/(F);
          end          
      end
end