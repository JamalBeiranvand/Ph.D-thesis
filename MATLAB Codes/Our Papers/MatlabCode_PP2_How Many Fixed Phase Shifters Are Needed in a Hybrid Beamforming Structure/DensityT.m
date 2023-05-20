clear
clc
N=17
if ~mod(ceil(N/2),2)
    n=1:ceil(N/4);
    Eta=2*sum(cos(pi/N.*(2*n-1)));
else    
    n=-floor(N/4):floor(N/4);
    Eta=sum(cos(2*pi*n./N));
end

S=pi*Eta^2;
NdistinctPoints(N)
D=NdistinctPoints(N)/S