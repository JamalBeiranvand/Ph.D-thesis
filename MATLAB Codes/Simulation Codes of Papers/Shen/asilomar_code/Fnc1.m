function V=Fnc1(obj,weight)
L = obj.numBS;
K = obj.numUser;
M = obj.numTxAnte;
N = obj.numRxAnte;
noise = obj.noise;
H = obj.chn;
maxPower = obj.maxPower;
association = obj.association;
schedule=zeros(L,K);
for i=1:L
    Im=find(association==i);
    BS(i)=length(Im);
    schedule(i,1:length(Im))=Im;
end
schedule;
% schedule = initScheduler(obj, weight);
numIter=100;
%initialize U & V
V = nan(M,K);
for i = 1:K
    V(:,i) = ones(M,1)*sqrt(maxPower(i)/M);
end
currentSINR=ZF(H,sum(maxPower),noise);
currentRate = sum(obj.bandwidth*log2(1+currentSINR))

% algorithm iterations
for iter = 1:numIter
    gamma = computeCurrentSINR1(obj, schedule, V);   
	Y= updateY( L, K, noise, H, schedule, V, weight, gamma);
    V= updateScheduleAndV(L, M, H, maxPower,schedule, weight, Y, gamma);
    norm(V,'fro')
    currentSINR=computeCurrentSINR1(obj, schedule, V);
    currentRate = sum(obj.bandwidth*log2(1+currentSINR))
end

end

%%
function [ currentSINR ] = computeCurrentSINR1(obj, schedule, V)

L = obj.numBS;
K = obj.numUser;
N = obj.numRxAnte;
H = obj.chn;
noise = obj.noise;
currentSINR = zeros(K,1);

for j = 1:L
    for s = 1:sum(schedule(j,:)~=0)
        i = schedule(j,s);
        A = norm(H(:,:,i,j)*V(:,i))^2;
        if A == 0
            continue
        end      
        B = noise;
        for n = 1:L
            for t = 1: sum(schedule(n,:)~=0)
                m = schedule(n,t);
                B = B + norm(H(:,:,i,j)*V(:,m))^2;
            end
        end
        
        currentSINR(i) = A/(B-A);
    end
end

end

%%
function [ Y ] = updateY( L, K, noise, H, schedule, V, weight, gamma)

Y = zeros(L,K);
        for i = 1:L
            for k=1:sum(schedule(i,:)~=0)
                m = schedule(i,k);
                B = noise;
                for j = 1:L
                    for t = 1:sum(schedule(j,:)~=0)
                        n = schedule(j,t);
                        B = B + H(:,:,m,j)*V(:,n)*V(:,n)'*H(:,:,m,j)';
                    end
                end
                A = sqrt(weight(m)*(1+gamma(m)))*H(:,:,m,i)*V(:,m);
                Y(i,m) = B \ A;
            end
        end
end
            
%%
function V  = updateScheduleAndV(L, M, H, maxPower,schedule, weight, Y, gamma)

        for i = 1:L
            for k=1:sum(schedule(i,:)~=0)
                m = schedule(i,k);
                B = zeros(M,M);
                for j = 1:L
                    for t = 1:sum(schedule(j,:)~=0)
                        n = schedule(j,t);
                        B = B + H(:,:,n,j)'*Y(j,n)*Y(j,n)'*H(:,:,n,j);
                    end
                end
                A = sqrt(weight(m)*(1+gamma(m)))*H(:,:,m,i)'*Y(i,m);
                tempZ = B \ A;
        
                if norm(tempZ)^2 <= maxPower(i)
                   V(:,m) = tempZ;
                else
                    % bisection search on opt mu
                    muRight = 1e-8;
                    while 1
                        tempZ = (B+muRight*eye(M))\A;
                        if norm(tempZ)^2 < maxPower(i)
                            V(:,m) = tempZ;
                            break 
                        end
                        muLeft = muRight;
                        muRight = muRight*5;
                    end
                end
            end
        end
end

function currentSINR=ZF(H,P,noise)
[~,M,K]=size(H);
Hzf=zeros(K,M);
for k=1:K
    Hzf(k,:)=H(1,:,k);
end
Fz=Hzf'*inv(Hzf*Hzf');
Fz=sqrt(P)*Fz/norm(Fz,'fro');
Pzf=norm(Fz,'fro')
diag(Hzf*Fz)
currentSINR=diag(Hzf*Fz)/noise;
end

% assignment problem