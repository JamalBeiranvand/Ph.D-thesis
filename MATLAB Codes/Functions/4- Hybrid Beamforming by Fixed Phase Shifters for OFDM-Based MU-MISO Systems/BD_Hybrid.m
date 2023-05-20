function  FBD= BD_Hybrid(Ht,WBB,WRF,FRF,P)
[NRF,Ns,K,F]=size(WBB);


for f=1:F
    %FBf=[];
%     for k=1:K
%         FBf=[FBf,FBB(:,:,k,f)];
%     end
    for k=1:K       
        W=WRF(:,:,k)*WBB(:,:,k,f);
        H(:,:,k,f)=W'*Ht(:,:,k,f)*FRF;
    end
end
[Nr,~,~,~]=size(H);

for f=1:F
    Hs=[];
    Ms=[];
    for k=1:K
        Hs=[Hs;H(:,:,k,f)];
    end
    Sig=[];
    for k=1:K
        Hj=H(:,:,k,f);
        Hhatj=Hs;
        Hhatj((k-1)*Nr+1:k*Nr,:)=[];
        Lhatj=rank(Hhatj);
        [~,~,Vhatj]=svd(Hhatj);
%         Vhatj1=Vhatj(:,1:Lhatj);
        Vhatj0=Vhatj(:,Lhatj+1:end);
        [~,S,Vj]=svd(Hj*Vhatj0);
        Sig=[Sig,diag(S(1:Ns,1:Ns)).'];
        Vj1=Vj(:,1:Ns);
        Mkf=[Vhatj0*Vj1];
        FBD(:,:,k,f)=sqrt(Ns*P)*Mkf/norm(FRF*Mkf,'fro');
        Ms=[Ms,Mkf];
    end    
%     G = func_water_filling_algorithm((1)*ones(K*Ns,1),Sig.^2,Ns*K*P);
%     MS=Ms*sqrt(G);
%     for k=1:K
%         FBD(:,:,k,f)=MS(:,(k-1)*Ns+1:k*Ns);
%     end
end
% X=reshape(Fopt,1,[]);

% Fopt=sqrt(K*Ns*F)*Fopt/norm(X,'fro');
% Fopt=sqrt(P)*Fopt/norm(X,'fro');


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