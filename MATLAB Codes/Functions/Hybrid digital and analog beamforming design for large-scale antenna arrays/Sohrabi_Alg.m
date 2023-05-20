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
function [WBB,WRF,FRF,FBB]=Sohrabi_Alg(H,NRF,Ns,P,Nq)
if nargin < 4
    error('At least 4 inputs are required');
end
if nargin == 4
    Nq='inf';
end
    [Mr,Mt]=size(H);
    FRF=ones(Mt,NRF);
    F1=H'*H;
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
    H_e = H*FRF/(Q12);
    [~,S,V] = svd(H_e);   % H = U*S*V'
    V1 = V(:,1:Ns);
    water_power = func_water_filling_algorithm((1)*ones(Ns,1),diag(S(1:Ns,1:Ns).^2),P);
    FBB = (Q12)\(V1*sqrt(water_power));
    
%    trace(Q*FBB*FBB')
    Vt=FRF*FBB;
    gamma2=1/Mr;
    F2=H*Vt*(Vt)'*H';
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
    J=(WRF)'*H*Vt*(Vt)'*(H)'*WRF+ (WRF)'*WRF;
    WBB=J\((WRF)'*H*Vt);
    
%     %%
%     function G=WaterFill(qkk,Pt)
%     [qkknew,indx]=sort(qkk,'descend');
%     N=length(qkknew);
%     t=1; %k is the first value for while
%     G=ones(N,1);% g is vector of gain at transmiter 
%     a=-1;% a is the first value for while
%     while a<0
%         Lambda= (N-t+1)/(Pt+sum(qkknew(1:(N-t+1))));
%         G(1:N-t+1)=1/Lambda-((qkknew(1:(N-t+1))));% comput the vector g
%         a=G(N-t+1);% compute the k th element of vector g
%         if a<0  % if the Kth element is negetive, it change to zero 
%             G(N-t+1)=0;
%             t=t+1; % compute the next element of vector g
%         end
%     end 
%     G(indx)=G;
%     end
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