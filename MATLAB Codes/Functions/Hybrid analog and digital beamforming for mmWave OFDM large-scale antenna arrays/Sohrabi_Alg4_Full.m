% Note: you need "Qtz" function to use this function

% This function is written according to
% "Hybrid beamforming and one-bit precoding for large-scale antenna arrays" 
% you can see the paper on this link: 
% https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=Hybrid+beamforming+and+one-bit+precoding+for+large-scale+antenna+arrays&btnG=
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
function [WBB,WRF,FRF,FBB]=Sohrabi_Alg4_Full(H,NRF,Ns,P,Nq)
if nargin < 4
    error('At least 4 inputs are required');
end
if nargin == 4
    Nq='inf';
end
    [Mr,Mt,F]=size(H);
    FRF=ones(Mt,NRF);
    F1=zeros(Mt);
    for f=1:F
        F1=F1+H(:,:,f)'*H(:,:,f);           
    end
    F1=F1/F;
        
    gamma2=P/(Mt*NRF);
    Delat=10;
    Rnew=0;
    while Delat>.1
        for j=1:NRF
            VRFj=FRF;
            VRFj(:,j)=[];
            Cj=eye(NRF-1)+gamma2*VRFj'*F1*VRFj;
            Gj=gamma2*F1-((gamma2^2)*F1*VRFj*(Cj\((VRFj)'*F1)));
            Gj=Gj - diag(diag(Gj)); 
            for i=1:Mt
                eta=Gj(i,:)*FRF(:,j);
                if eta~=0
                    if nargin == 5
                        FRF(i,j)=Qtz(eta/abs(eta),Nq);
                    else
                        FRF(i,j)=eta/abs(eta);
                    end
                else
                    FRF(i,j)=1;
                end
            end
        end
        Rold=Rnew;
        Rnew=log2(det(eye(NRF)+gamma2*FRF'*F1*FRF));
        Delat=Rnew-Rold;
    end
    Q = (FRF'*FRF);
    Q12 = sqrtm(Q);
    
    for f=1:F
        H_e = H(:,:,f)*FRF;
        [~,S,V] = svd(H_e/Q12);   % H = U*S*V'
        V1 = V(:,1:Ns);
        water_power = func_water_filling_algorithm((1)*ones(Ns,1),diag(S(1:Ns,1:Ns).^2),P);
        FBB(:,:,f) = (Q12)\(V1*sqrt(water_power));
    end
        
    for f=1:F
        Vt=FRF*FBB(:,:,f);
        Fhat(:,:,f)=H(:,:,f)*Vt*(Vt)'*H(:,:,f)';           
    end
    
    F2=zeros(Mr);
    for f=1:F
        F2=F2+Fhat(:,:,f);           
    end
    F2=F2/F;
%    trace(Q*FBB*FBB')
%     Vt=FRF*FBB;
    gamma2=1/Mr;
    WRF=ones(Mr,NRF);
    Delat=10;
    Rnew=0;
    while Delat>.1
        for j=1:NRF
            WRFj=WRF;
            WRFj(:,j)=[];
            Cj=eye(NRF-1)+gamma2*WRFj'*F2*WRFj;
            Gj=gamma2*F2-(gamma2^2)*F2*WRFj*(Cj\((WRFj)'*F2));
            Gj=Gj - diag(diag(Gj)); 
            for i=1:Mr
                eta=Gj(i,:)*WRF(:,j);
                if eta~=0
                   if nargin == 5
                        WRF(i,j)=Qtz(eta/abs(eta),Nq);
                   else
                        WRF(i,j)=eta/abs(eta);
                   end
                else
                    WRF(i,j)=1;
                end
            end
        end
        Rold=Rnew;
        Rnew=log2(det(eye(NRF)+gamma2*WRF'*F2*WRF));
        Delat=Rnew-Rold;
    end
    
    WBB=[];
    for f=1:F
        Vt=FRF*FBB(:,:,f);
        J=(WRF)'*H(:,:,f)*Vt*(Vt)'*(H(:,:,f))'*WRF+ (WRF)'*WRF;
        WBB(:,:,f)=J\((WRF)'*H(:,:,f)*Vt);
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function water_power = func_water_filling_algorithm(SNR,lambda,Trans_Power)

        Noise_Power= (1./(SNR.*lambda));
        Number_Channel= length(Noise_Power) ;
        [S_Number, dt]=sort(Noise_Power);
        sum(Noise_Power);
        for p=length(S_Number):-1:1
            T_P=(Trans_Power+sum(S_Number(1:p)))/p;
            Input_Power=T_P-S_Number;
            Pt=Input_Power(1:p);
            if(Pt(:)>=0)
                break
            end
        end
        Allocated_Power=zeros(1,Number_Channel);
        Allocated_Power(dt(1:p))=Pt;
        water_power = diag(Allocated_Power);
        end
end