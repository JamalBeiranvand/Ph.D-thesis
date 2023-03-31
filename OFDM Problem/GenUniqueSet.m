function [Setb,Sb,Set,S]=GenUniqueSet(N)
P=unique(factor(N));
Xi=N-N/min(P);
SET=1:N;
Validpoint=cell(Xi,1);
for g=1:Xi
    Selection=nchoosek(SET,g);   
    Validpoint{g,1}=Selection;
end

S=false(1,N);
t=2;
for i=1:Xi
    L=size(Validpoint{i,1},1);
    for j=1:N
        S(t:t+L-1,j)=logical(sum(Validpoint{i,1}==j,2));        
    end
    t=t+L;
end

for j=1:length(P)
    for i=1:N/P(j)
        A=i:N/P(j):N;
        Sa=5*ones(1,N);
        Sa(A)=1;
        S(sum(S==Sa,2)==length(A),:)=[];
    end
end

if ~mod(N,2)&& length(P)>1
    P(P==2)=[];
    for j=1:length(P)
        for i=1:N/P(j)
            A=i:N/P(j):N;
            for k=((P(j)+1)/2):(P(j)-1)
                T=nchoosek(A,k);
                for z=1:size(T,1)
                    Sa=5*ones(1,N);
%                     Sa(A)=0;
                    Sa(T(z,:))=1;
                    Cond1=(sum(S==Sa,2)==k);
                    Cond2=(sum(S==Sa,2)==sum(S,2));
                    S(Cond1&Cond2,:)=[];
                end
            end
        end
    end
end

if N==12
    for j=1:length(P)
        for i=1:N/P(j)
            A=i:N/P(j):N;
            if i<3
                 B=i+2:N/P(j):N;
            else
                 B=i-2:N/P(j):N;
            end
            for k=((P(j)+1)/2):(P(j)-1)
                T=nchoosek(A,k);
                for z=1:size(T,1)
                    Sa=5*ones(1,N);
                    Sa(A)=0;
                    Sa(T(z,:))=1;
                    Sa(B)=0;
                    Cond1=(sum(S==Sa,2)==6);
                    S(Cond1,:)=[];
                end
            end
        end
    end
end

phases=exp(1j*2*pi*(0:N-1)/N);
Set=zeros(size(S,1),1);
for i=2:size(S,1)
 Set(i)=sum(phases(S(i,:)));
end

indx=angle(Set)>=-10e-6 & angle(Set)<2*pi/N-10e-6;
Setb=Set(indx);       % points in basis set 
Sb=S(indx,:); 