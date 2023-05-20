%% greedy Algorithm
function [ComSet,Hc]=ComAntSelection_Greedy(H,NCom)
  [Nt,K]=size(H);
    if Nt<K
        error('size of the channel matrix is not true')
    end
    
    X=H;
    ComSet=1:Nt;
    for i=1:Nt-NCom
        D=[];
        for n=1:length(ComSet)
            hh=H(ComSet(n),:).'*conj(H(ComSet(n),:));
            D(n)=trace(inv(X.'*conj(X)-hh));
        end 
        N=D==min(D);
        X(N,:)=[];
        ComSet(N)=[];
    end 
    Hc=H(ComSet,:);
    