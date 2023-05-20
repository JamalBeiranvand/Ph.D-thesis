function  [Wopt,FBD]= BD(H,Ns,P)
% [Nr,Nt,K,F]=size(Ht);
% 
%  
%    for f=1:F
%        Ff=[];
%        for k=1:K
%            [Wop,~,Fop]=svd(Ht(:,:,k,f));
%             Wopt(:,:,k,f)=Wop(:,1:Ns);
%             Ff=[Ff,Fop(:,1:Ns)];
%        end
%        [UF,~,~]=svd(Ff);
%        Hf(:,:,f)=UF(:,1:K*Ns);
%    end
%    
%    for f=1:F
%         for k=1:K       
%             H(:,:,k,f)=Wopt(:,:,k,f)'*Ht(:,:,k,f)*Hf(:,:,f);
%         end
%    end
[Nr,Nt,K,F]=size(H);

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
        Ms=[Ms,Mkf];
    end  
    G = func_water_filling_algorithm((1)*ones(K*Ns,1),Sig,Ns*K*P);
    MS=Ms*sqrt(G);
    for k=1:K
        FBD(:,:,k,f)=MS(:,(k-1)*Ns+1:k*Ns);
    end
end

for f=1:F
    for k=1:K
        [W,~,~]=svd(H(:,:,k,f)*FBD(:,:,k,f));
        Wopt(:,:,k,f)=W(:,1:Ns);
    end
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