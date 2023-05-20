function [KSets,Nuk]=NonSharedFixed(Nu,K)
   Nuk=floor(Nu/K)*ones(K,1);
   k=1;
   while sum(Nuk)<Nu
        Nuk(k)=Nuk(k)+1;
        k=k+1;
   end
   KSets=zeros(max(Nuk),K);
   nt=1;
   for k=1:K
       for nuk=1:Nuk(k)
           KSets(nuk,k)=nt;
           nt=nt+1;
       end
   end
end
