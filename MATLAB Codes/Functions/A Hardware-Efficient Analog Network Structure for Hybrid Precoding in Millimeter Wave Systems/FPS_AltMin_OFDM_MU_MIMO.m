function [WBB,WRF,FRF,FBB] = FPS_AltMin_OFDM_MU_MIMO(FBD,W,NrfTx,NrfRx,P,Np)

  [Nt,Ns,K,F]=size(FBD);
   Fopt=[];  
   for kk=1:K
        for f=1:F
            Fopt=[Fopt,FBD(:,:,kk,f)];
        end 
   end
    
C=Phases(Np,NrfTx)./sqrt(Np);

[~,~,V] = svd(Fopt);
FD = V(:,1:NrfTx)';
indx=0;
mynorm = [Inf,0];
% while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
while indx<10
   [alpha, value, S] = alpha_opt_new(real(Fopt*FD'*C'));
   
    S = reshape(S,[Nt,Np*NrfTx]);
    mynorm(1) = value + norm(imag(Fopt*FD'*C'),'fro')^2;

    [U,~,V] = svds(alpha*Fopt'*S*C,NrfTx);
    FD = V*U';
    mynorm(2) = norm(Fopt*FD'*C' - alpha*S,'fro')^2;
    indx=1+indx;
end
FD = alpha*FD;
FRF = S*C;
for kk=1:K
    for f=1:F
        St=(kk-1)*F+(f-1)*Ns+1;
        Ed=(kk-1)*F+f*Ns;
        FBB(:,:,kk,f)=FD(:,St:Ed);
        FBB(:,:,kk,f) = sqrt(Ns*P) * FBB(:,:,kk,f) / norm(FRF * FBB(:,:,kk,f),'fro');
    end
end

%% Reciver side
C=Phases(Np,NrfRx)./sqrt(Np);
[Nr,~,~,~]=size(W);
for kk=1:K
    Wopt=[];
    for f=1:F
        Wopt=[Wopt,W(:,:,kk,f)];
    end 
    [~,~,V] = svd(Wopt);
    WD = V(:,1:NrfRx)';
    indx=0;
    mynorm = [Inf,0];
    % while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
    while indx<10
       [alpha, value, S] = alpha_opt_new(real(Wopt*WD'*C'));

        S = reshape(S,[Nr,Np*NrfRx]);
        mynorm(1) = value + norm(imag(Wopt*WD'*C'),'fro')^2;

        [U,~,V] = svds(alpha*Wopt'*S*C,NrfRx);
        WD = V*U';
        mynorm(2) = norm(Wopt*WD'*C' - alpha*S,'fro')^2;
        indx=1+indx;
    end
    WD = alpha*WD;
    WRF(:,:,kk) = S*C;
    for f=1:F
        WBB(:,:,kk,f)=WD(:,(f-1)*Ns+1:f*Ns);
    end
end

%% alpha_opt_new 
    function [Alpha,objt,Sss]=alpha_opt_new(x)
    [x,ord] = sort([x(:);0],'ascend');
    n = length(x);
    ave = zeros(n-1,1);
    ave(1) = x(1);
    ave(end) = x(end);
    k = find(x==0,1);
    for i = 2:k-1
        ave(i) = (ave(i-1)*(i-1)+x(i))/i;
    end
    for i = n-2:-1:k
        ave(i) = (ave(i+1)*(n-i-1)+x(i+1))/(n-i);
    end

    idx = find(ave>=2*x(1:n-1) & ave<=2*x(2:n));

    obj = zeros(2,length(idx));
    for i = 1:length(idx)
        Alpha = ave(idx(i));
        if(Alpha<=0)
            s(:,i) = double(x-Alpha/2<=0);
            obj(:,i) = [Alpha;norm(x-Alpha*s(:,i),2)^2];
        else
            s(:,i) = double(x-Alpha/2>0);
            obj(:,i) = [Alpha;norm(x-Alpha*s(:,i),2)^2];
        end
    end
    [objt,loc] = min(obj(2,:));
    % objt = objt-norm(x,2)^2;
    Alpha = obj(1,loc);
    temp = s(:,loc);
    Sss(ord) = temp;
    Sss(end)=[];
    end

end