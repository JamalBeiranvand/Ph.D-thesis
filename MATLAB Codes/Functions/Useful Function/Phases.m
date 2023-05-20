function [C,c]=Phases(N,Ns) 
c=exp(1j*2*pi*(0:N-1)/N).';
C=[];
for i=1:Ns
    C=blkdiag(C,c);
end