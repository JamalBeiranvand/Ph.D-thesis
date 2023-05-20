 % Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
 % google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
 % Date : 29/05/2020
function R=SumRate(H,F,G,varargin)
%% Check Inputs
    p = inputParser; 
    DefaultNoiseVar=1;
    DefaultLink='Down';
    validScalarPosNum= @ (x) isnumeric(x) && isscalar(x) && (x > 0);
    addRequired(p,'H');
    addRequired(p,'F');
    addRequired(p,'G');
    addParameter(p,'NoiseVar',DefaultNoiseVar,validScalarPosNum);
    
    addParameter(p,'Link',DefaultLink);
    parse(p,H,F,G,varargin{:});
    Sigma2=p.Results.NoiseVar;
    Link=p.Results.Link;
%% Main Code
if strcmp(Link,'Down')
    [K,~]=size(H);
    A=H*F*G;
    Rk=zeros(K,1);
    for k=1:K
        Num=abs(A(k,k))^2;
        Dem=Sigma2+sum(abs(A(:,k)).^2)-Num;
        Rk(k)=log2(1+Num/Dem);
    end
    R=sum(Rk);
end

if strcmp(Link,'Up')
    W=F;
    [~,K]=size(H);
    A=W'*H*G;
    Rk=zeros(K,1);
    for k=1:K
        Num=abs(A(k,k))^2;
        Wk=sum(abs(W(:,k)).^2);
        Dem=Wk+sum(abs(A(:,k)).^2)-Num;
        Rk(k)=log2(1+Num/Dem);
    end
    R=sum(Rk);
end
end