function [Setb,Sb,Set,S]=GenSet(N)
SET=1:N;
Validpoint=cell(N/2,1);
for g=1:N/2
    Selection=nchoosek(SET,g);
    C=1;
    i=1;
    t=1;
    while C>.1    
        sub=Selection(i,:);
            for j=1:N/2
                if (sum(sub==j)&& sum(sub==j+N/2))
                  Selection(i,:)=[];
                  i=i-1;
                  break
                end    
            end
        if i==length(Selection)
            C=0;
        end 
        i=i+1;
        t=1;
    end
    Validpoint{g,1}=Selection;
end
phases=exp(1j*2*pi*(0:N-1)/N).';
Set=[0;phases];
 S=false(1,N);
t=2;
for i=1:N/2
    L=length(Validpoint{i,1});
    for j=1:N
        S(t:t+L-1,j)=logical(sum(Validpoint{i,1}==j,2));        
    end
    t=t+L;
end

for i=2:N/2
  Snew=phases(Validpoint{i,1});
  Set=[Set;sum(Snew,2)];
end 
% 
% extract basis set
    indx=angle(Set)>=0 & angle(Set)<=2*pi/N;
    Setb=Set(indx);       % points in basis set 
    Ssb=S(indx,:);       % switchs states of points in basis set
% sort basis set
    [~,idx]=sort(abs(Setb));
    Setb=Setb(idx);      % sorted points
    Sb=Ssb(idx,:);% order swithes according sorted points
% Generate C matrix 
