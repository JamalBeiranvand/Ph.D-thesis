%Z = Laprand(Mu,b,M,N) generates a M by N array containing random numbers
%from the Laplacian distribution with  parameters mu  and  b
function Z=Laprand(Mu,b,M,N)
if nargin == 2
    M=1;
    N=1;
end
if nargin == 3
    N=M;
end
if nargin < 2
    error('At least 2 inputs are required');
end
X1=exprnd(b,M,N);
X2=exprnd(b,M,N);
Z=Mu+X1-X2;
