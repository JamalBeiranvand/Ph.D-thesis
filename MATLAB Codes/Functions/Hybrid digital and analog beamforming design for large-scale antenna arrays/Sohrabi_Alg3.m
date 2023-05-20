function [VRF,VD]= Sohrabi_Alg3(H,NRF,Ptotal,Nq)
% H is a K*Mt chanel matrix
% NRF is Number of RF chain at BS
% Pt is the total power budget at the BS
% power of noise Sigma=1;
if nargin == 3
    Nq='inf';
end
if nargin == 4
    set=2*pi*(0:Nq-1)/Nq;
end
[K,Mt]=size(H);
% 1.Start with a feasible VRF and P = IK.
    g=ones(K,1);
    P=diag(g);                              %% power matrix
    VRF=exp(1j*2*pi*randn(Mt,NRF));                       %% initial VRF  
    for Iter=1:2     %% 10 iteration for coveragence of the overal algotithm
        for It=1:5     %% 10 iteration for coveragence of VRF
            % 2. j=1:NRF do
            for j=1:NRF
                % 3. Calculate Aj
                VRFj=VRF;
                VRFj(:,j)=[];                 %% delet j-th colummn
                g(g>0)=g(g>0).^(-1/2);        %% because of 0^(-1/2)=NAN;
                P12=diag(g);                  %% P12=P^(-1/2)
                Aj=P12*H*VRFj*(VRFj')*H'*P12; %% acording line 3 of Algorithm 3
                % compute Bj and Dj 
                Hhat=P12*H;                   %% from P.39 under equation (2.32)
                Aj2=inv(Aj)*inv(Aj);          %% Aj2=Aj^(-2);
                Bj=Hhat'*Aj2*Hhat;            %% from P.40 under equation (2.35)
                Dj=Hhat'*(Aj\Hhat);      %% from P.40 under equation (2.35)
                
                % 4. for i=1:Mt do
                for i=1:Mt
                    B=Bj;                     
                    B(i,:)=0;
                    B(:,i)=0;
                    D=Dj;
                    D(i,:)=0;
                    D(:,i)=0;

                    ZetaBij=Bj(i,i)+2*real(VRF(:,j)'*B*VRF(:,j)); %% from (2.35a)
                    ZetaDij=Dj(i,i)+2*real(VRF(:,j)'*D*VRF(:,j)); %% from (2.35b)
                    EtaBij =Bj(i,:)*VRF(:,j)-Bj(i,i)*VRF(i,j);     %% from (2.35c)
                    EtaDij =Dj(i,:)*VRF(:,j)-Dj(i,i)*VRF(i,j);     %% from (2.35d)

                    cij=(1+ZetaDij)*EtaBij-ZetaBij*EtaDij;        %% from P.42 under equation (2.42b) 
                    zij=imag(2*conj(EtaBij)*EtaDij);              %% from P.42 under equation (2.42b)
                    
                    if strcmp(Nq,'inf')
                        % compute Phi from (2.43)
                        if real(cij)>=0             
                            Phiij=asin(imag(cij)/abs(cij));
                        else
                            Phiij=pi-asin(imag(cij)/abs(cij));
                        end

                        % compute Theta 1 and 2 
                        Theta1ij=  -Phiij+asin(zij/abs(cij));   %% from (2.42a)
                        Theta2ij=pi-Phiij-asin(zij/abs(cij));   %% from (2.42b)

                        % compute f(VRF(i,j)) for theta 1 and 2
                        f1=Fhat(Aj,ZetaBij,ZetaDij,EtaBij,EtaDij,Theta1ij);
                        f2=Fhat(Aj,ZetaBij,ZetaDij,EtaBij,EtaDij,Theta2ij);

                        % find Theta Optimum
                        if abs(f1)<abs(f2)
                            ThetaijOpt=Theta1ij;
                        else
                            ThetaijOpt=Theta2ij;
                        end
                        % 8: Set VRF(i,j)=e^(-j*Theta)
                        VRF(i,j)=exp(-1j*ThetaijOpt);
                    else
                        for s=1:length(set)
                            THeta=set(s);
                            fhat(s)=Fhat(Aj,ZetaBij,ZetaDij,EtaBij,EtaDij,THeta);
                        end
                        Nopt=find(abs(fhat)==min(abs(fhat)));
                        ThetaijOpt=set(Nopt(1));
                        VRF(i,j)=exp(-1j*ThetaijOpt);
                    end

                end  % 9: end for i=1:Mt     
            end      % 10: end 10 end for j=1:NRF
        end          % 11: end 10 iteration for convergence analog precoder
        
        % 12: Find P=diag[p1,...,PK] using waterfilling as in 2.29
        VDhat=VRF'*H'/(H*VRF*(VRF')*H');  %% from page 38
        qkk=diag((VDhat')*(VRF)'*VRF*VDhat);
        g=WaterFill(qkk,Ptotal); %% the WaterFill function is in the end of the algorithm
        P=diag(g./qkk);
        
    end % 13: end 10 iteration for convergence of the overal algotithm
    
    % 14. Set VD
    VD=VRF'*H'*((H*VRF*(VRF')*H')\P^(1/2)); 

% %% waterfilling function
    function G=WaterFill(qkk,Pt)
    [qkknew,indx]=sort(qkk,'descend');
    N=length(qkknew);
    t=1; %k is the first value for while
    G=ones(N,1);% g is vector of gain at transmiter 
    a=-1;% a is the first value for while
    while a<0
        Lambda= (N-t+1)/(Pt+sum(qkknew(1:(N-t+1))));
        G(1:N-t+1)=1/Lambda-((qkknew(1:(N-t+1))));% comput the vector g
        a=G(N-t+1);% compute the k th element of vector g
        if a<0  % if the Kth element is negetive, it change to zero 
            G(N-t+1)=0;
            t=t+1; % compute the next element of vector g
        end
    end 
    G(indx)=G;
    end
    %Fhat
    function F=Fhat(Aj,ZetaBij,ZetaDij,EtaBij,EtaDij,Theta)
        F=trace(inv(Aj))-((ZetaBij+2*real(conj(exp(-1j*Theta))*EtaBij))/(1+ZetaDij+2*real(conj(exp(-1j*Theta))*EtaDij)));
    end

% function water_power = func_water_filling_algorithm(SNR,lambda,Trans_Power)
% 
%         Noise_Power= (1./(SNR.*lambda));
%         Number_Channel= length(Noise_Power) ;
%         [S_Number, dt]=sort(Noise_Power);
%         sum(Noise_Power);
%         for p=length(S_Number):-1:1
%             T_P=(Trans_Power+sum(S_Number(1:p)))/p;
%             Input_Power=T_P-S_Number;
%             Pt=Input_Power(1:p);
%             if(Pt(:)>=0)
%                 break
%             end
%         end
%         Allocated_Power=zeros(1,Number_Channel);
%         Allocated_Power(dt(1:p))=Pt;
%         water_power = diag(Allocated_Power);
%         end
end
