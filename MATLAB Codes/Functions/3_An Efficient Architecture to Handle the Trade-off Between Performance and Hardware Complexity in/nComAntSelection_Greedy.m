function [KSets,d]=nComAntSelection_Greedy(H,ComSet)
    [Nt,K]=size(H);
    if Nt<K
        error('size of the channel matrix is not true')
    end
    
    
    if size(ComSet,1)<size(ComSet,2)
        ComSet=ComSet.';
    end
    
    Hc=H(ComSet,:);
    NonSharedSet=1:Nt;
    NonSharedSet(ComSet)=[];
    NC=length(ComSet);
    Nu=Nt-NC;
   
   Nuk=floor(Nu/K)*ones(K,1).';
   k=1;
   while sum(Nuk)<Nu
        Nuk(k)=Nuk(k)+1;
        k=k+1;
   end
   
    SetK=[repmat(ComSet,1,K);zeros(max(Nuk),K)];
    sk=NC*ones(1,K);
    d=diag(inv(Hc.'*conj(Hc))).';
     D=[];
    for n=1:length(NonSharedSet)
        for k=1:K
            Hk=[H(NonSharedSet(n),:);H(SetK(1:sk(k),k),:)];
            IvMx=inv(Hk.'*conj(Hk));
            D(n,k)=IvMx(k,k);                
        end
    end
        
    for i=1:Nu-1     
        DeltaD=repmat(d,length(NonSharedSet),1)-D;       
        [DeltaDS,IndxDeltaDS]=sort(DeltaD,'descend');        
        [~,Ks]=sort(DeltaDS(1,:),'descend');
        % Constraint on Nuk        
        KK=1:K;
        KK(sk~=(Nuk+NC))=[];
        for j=1:length(KK)
            Ks(Ks==KK(j))=[];
        end
        % find user
        k=Ks(1);
        IndxDeltaDS(1,k);
        Ns=NonSharedSet(IndxDeltaDS(1,k));
        %  add a antenna to k-th set
        sk(k)=sk(k)+1;
        SetK(sk(k),k)=Ns;
        % update Determinants
        Hk=H(SetK(1:sk(k),k),:);
        IvMx=inv(Hk.'*conj(Hk));
        d(k)=IvMx(k,k); 
        
         % delete allocted antenna      
        IndxNa=find(NonSharedSet==Ns);
        NonSharedSet(IndxNa)=[];
        D(:,k)=0;
        D(IndxNa,:)=[];        
        for n=1:length(NonSharedSet)
            Hk=[H(NonSharedSet(n),:);H(SetK(1:sk(k),k),:)];
            IvMx=inv(Hk.'*conj(Hk));
            D(n,k)=IvMx(k,k);                
        end     
    end
    k=find(sk~=(Nuk+NC));
    sk(k)=sk(k)+1;
    SetK(sk(k),k)=NonSharedSet;
    KSets=SetK(NC+1:end,:);
end
