function [WRF,WBB,FRF,FBB] = FPS_AltMin(H,NRF,Ns,Np)
[U,S,V] = svd(H);
Fopt = V(:,1:Ns);  
Wopt = U(:,1:Ns);
C=Phases(Np,NRF)./sqrt(Np);

[Nt, Ns] = size(Fopt);
[~,~,V] = svd(Fopt);
FBB = [V';zeros(NRF-Ns,Ns)];
indx=0;
mynorm = [Inf,0];
% while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
while indx<10
   [alpha, value, S] = alpha_opt_new(real(Fopt*FBB'*C'));
   
    S = reshape(S,[Nt,Np*NRF]);
    mynorm(1) = value + norm(imag(Fopt*FBB'*C'),'fro')^2;

    [U,~,V] = svds(alpha*Fopt'*S*C,NRF);
    FBB = V*U';
    mynorm(2) = norm(Fopt*FBB'*C' - alpha*S,'fro')^2;
    indx=1+indx;
end
FBB = alpha*FBB;
FRF = S*C;
FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');

%% Reciver side
[Nr, Ns] = size(Wopt);
NRF = size(C,2);
Np = size(C,1)/NRF;
[~,~,V] = svd(Wopt);
WBB = [V';zeros(NRF-Ns,Ns)];
indx=0;
mynorm = [Inf,0];
% while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
while indx<10
   [alpha, value, S] = alpha_opt_new(real(Wopt*WBB'*C'));
   
    S = reshape(S,[Nr,Np*NRF]);
    mynorm(1) = value + norm(imag(Wopt*WBB'*C'),'fro')^2;

    [U,~,V] = svds(alpha*Wopt'*S*C,NRF);
    WBB = V*U';
    mynorm(2) = norm(Wopt*WBB'*C' - alpha*S,'fro')^2;
    indx=1+indx;
end
WBB = alpha*WBB;
WRF = S*C;
WBB = sqrt(Ns) * WBB / norm(WRF * WBB,'fro');


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
