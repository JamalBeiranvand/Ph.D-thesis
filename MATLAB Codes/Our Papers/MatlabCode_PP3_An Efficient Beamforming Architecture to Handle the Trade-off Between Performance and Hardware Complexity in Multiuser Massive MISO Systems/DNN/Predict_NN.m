function [ComSet,KSets,NN_pre]=Predict_NN(NN_Probability,Ncom,NAntenna_User)
    [Nt,K]=size(NN_Probability);
    ProbabilitySum=sum(NN_Probability,2);
    [~,Indx]=sort(ProbabilitySum);
    NN_pre=zeros(size(NN_Probability));
    KSets=zeros(max(NAntenna_User),K);
    Nuk=zeros(K,1);
    for n=1:Nt-Ncom
        [~,k]=max(NN_Probability(Indx(n),:));
        NN_pre(Indx(n),k)=1;
        Nuk(k)=Nuk(k)+1;
        KSets(Nuk(k),k)=Indx(n);            
        if Nuk(k)==NAntenna_User(k)
            NN_Probability(:,k)=-1;
        end
    end
    NN_pre(Indx(end-Ncom+1:end),:)=1;
    ComSet=Indx(end-Ncom+1:end);
end