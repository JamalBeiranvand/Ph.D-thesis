
function G=WaterFilling(Lambda,Pt)
[lambda,indx]=sort(Lambda,'descend');
N=length(lambda);
p=1; %k is the first value for while
G=ones(N,1);% g is vector of gain at transmiter 
a=-1;% a is the first value for while
while a<0
    mu= (N-p+1)/(Pt+sum(1./lambda(1:(N-p+1))));
    G(1:N-p+1)=1/mu-(1./(lambda(1:(N-p+1))));% comput the vector g
    a=G(N-p+1);% compute the k th element of vector g
    if a<0  % if the Kth element is negetive, it change to zero 
        G(N-p+1)=0;
        p=p+1; % compute the next element of vector g
    end
end 
G(indx)=G;