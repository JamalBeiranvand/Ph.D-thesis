clear
clc
K=2;
Nt=21;
SNR=0;% in dB
Snr=10^(SNR/10);
%% 
H=sqrt(1/2)*(randn(K,Nt)+1j*randn(K,Nt)).';
[Hsorted,INDX]=sort(abs(H));
Antenna=[];
for i=0:Nt/2
    ComAntenna=INDX(end-i,2);
    N=find(INDX(:,1)==ComAntenna);
    if N>Nt/2
        Antenna=[Antenna ComAntenna];
    end
end

for S=1:length(Antenna)
    Index=INDX;
    ComAntenna=Antenna(S)
    Index(Index==ComAntenna)=[];
    Index=reshape(Index,[],K);
    Frf=zeros(Nt,K);
    k=1;
    for n=1:Nt-1
        Frf(Index(end,k),k)=1;
        Index(Index==Index(end,k))=[];
        Index=reshape(Index,[],K);
        if k==1
            k=2;
        else
            k=1;
        end
    end
    Frf(Frf>0)=conj(H(Frf>0));
    FRF=Frf;
    
    IntF=H.'*Frf;
    Frf(ComAntenna,1)=-IntF(2,1)/H(ComAntenna,2);
    Frf(ComAntenna,2)=-IntF(1,2)/H(ComAntenna,1);
    PFrf=sum(sum(abs(Frf).^2));
    Frf=sqrt(Snr)*Frf/sqrt(PFrf);
    SINR=abs(H.'*Frf)
    
    %% water
    Lambda=H.*FRF; 
    lambda=Lambda;
    L1=Lambda(Lambda(:,1)>0,1);
    F1=WaterFilling(L1,Snr);
    L2=Lambda(Lambda(:,2)>0,2);
    F2=WaterFilling(L2,Snr); 
    
    Lambda(Lambda(:,1)>0,1)=F1;
    Lambda(Lambda(:,2)>0,2)=F2;
    
    FRFn=FRF.*Lambda;
    Int=H.'*FRFn;
    F1=WaterFilling(L1,Snr-Int(1,2));        
    F2=WaterFilling(L2,Snr-Int(2,1));  
    lambda(lambda(:,1)>0,1)=F1;
    lambda(lambda(:,2)>0,2)=F2;
    FRF=FRF.*lambda;
    
    IntF=H.'*FRF;
    FRF(ComAntenna,1)=-IntF(2,1)/H(ComAntenna,2);
    FRF(ComAntenna,2)=-IntF(1,2)/H(ComAntenna,1);
    PFrf=sum(sum(abs(FRF).^2));
    FRF=sqrt(Snr)*FRF/sqrt(PFrf);
    
    SINR2=abs(H.'*FRF)
end
%% ZF precoder
H=H.';
Azf=H'*inv(H*H');
PZf=sum(sum(abs(Azf).^2));
Azf=sqrt(Snr)*Azf/sqrt(PZf); 
ZFSINR=abs(H*Azf)
PZf=sum(sum(abs(Azf).^2));
    
%%MMSE Precoder
Ammse=H'/(H*H'+K/Snr*eye(K));
Pmmse=sum(sum(abs(Ammse).^2));
Ammse=sqrt(Snr)*Ammse/sqrt(Pmmse);
MMSESINR=abs(H*Ammse)