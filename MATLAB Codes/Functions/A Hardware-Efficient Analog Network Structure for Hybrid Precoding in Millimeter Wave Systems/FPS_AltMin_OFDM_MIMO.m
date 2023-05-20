function [WBB,WRF,FRF,FBB] = FPS_AltMin_OFDM_MIMO(H,NRF,Ns,P,Np)
[Nr,Nt,F]=size(H);
   Fopt=[];  
   Wopt=[];
    for f=1:F
        [Wop,S,Fop]=svd(H(:,:,f));
        Fopt=[Fopt,Fop(:,1:Ns)];
        Wopt=[Wopt,Wop(:,1:Ns)];
    end   
    
C=Phases(Np,NRF)./sqrt(Np);

[~,~,V] = svd(Fopt);
FD = V(:,1:NRF)';
indx=0;
mynorm = [Inf,0];
% while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
while indx<10
   [alpha, value, S] = alpha_opt_new(real(Fopt*FD'*C'));
   
    S = reshape(S,[Nt,Np*NRF]);
    mynorm(1) = value + norm(imag(Fopt*FD'*C'),'fro')^2;

    [U,~,V] = svds(alpha*Fopt'*S*C,NRF);
    FD = V*U';
    mynorm(2) = norm(Fopt*FD'*C' - alpha*S,'fro')^2;
    indx=1+indx;
end
FD = alpha*FD;
FRF = S*C;
for f=1:F
 FBB(:,:,f)=FD(:,(f-1)*Ns+1:f*Ns);
 FBB(:,:,f) = sqrt(P) * FBB(:,:,f) / norm(FRF * FBB(:,:,f),'fro');
end

%% Reciver side
[~,~,V] = svd(Wopt);
WD = V(:,1:NRF)';
indx=0;
mynorm = [Inf,0];
% while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
while indx<10
   [alpha, value, S] = alpha_opt_new(real(Wopt*WD'*C'));
   
    S = reshape(S,[Nr,Np*NRF]);
    mynorm(1) = value + norm(imag(Wopt*WD'*C'),'fro')^2;

    [U,~,V] = svds(alpha*Wopt'*S*C,NRF);
    WD = V*U';
    mynorm(2) = norm(Wopt*WD'*C' - alpha*S,'fro')^2;
    indx=1+indx;
end
WD = alpha*WD;
WRF = S*C;
for f=1:F
 WBB(:,:,f)=WD(:,(f-1)*Ns+1:f*Ns);
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