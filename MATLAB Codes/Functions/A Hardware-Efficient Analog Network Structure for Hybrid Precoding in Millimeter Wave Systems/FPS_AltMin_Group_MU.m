function [FRF,FBB] = FPS_AltMin_Group_MU(H,NRF,K,Np,Ng)
Fopt = ZF_Down(H); 
C=Phases(Np,NRF)./sqrt(Np);
Chat=Phases(Np,NRF/Ng)./sqrt(Np);
[Nt, K] = size(Fopt);
[~,~,V] = svd(Fopt);
FBB = [V';zeros(NRF-K,K)];
Sm=[];
FB=[];
for g=1:Ng
    Fi=Fopt((g-1)*(Nt/Ng)+1:g*(Nt/Ng),:);
    Bi=FBB((g-1)*(NRF/Ng)+1:g*(NRF/Ng),:);
        indx=0;
        mynorm = [Inf,0];
        % while (isempty(mynorm) || abs( mynorm(1) - mynorm(2) ) > 1e-5)
        while indx<10
           [alpha, value, S] = alpha_opt_new(real(Fi*Bi'*Chat'));

            S = reshape(S,[Nt/Ng,Np*NRF/Ng]);
            mynorm(1) = value + norm(imag(Fi*Bi'*Chat'),'fro')^2;

            [U,~,V] = svds(alpha*Fi'*S*Chat,NRF/Ng);
            Bi = V*U';
            mynorm(2) = norm(Fi*Bi'*Chat' - alpha*S,'fro')^2;
            indx=1+indx;
        end
        Bi = alpha*Bi;        
        Sm=blkdiag(Sm,S);
        FB=[FB;Bi];
end

FRF = Sm*C;
FBb = sqrt(K) * FB / norm(FRF * FB,'fro');


He= H*FRF*FBb; 
FBD=He'*inv(He*He');
FBB=FBb*FBD/norm(FRF * FBb*FBD,'fro');


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
