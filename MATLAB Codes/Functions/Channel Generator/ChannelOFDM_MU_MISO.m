function [H,At,Ar,alpha]=ChannelOFDM_MU_MISO(Tx,K,Ncl,F,varargin)
    
% Check Inputs
    p = inputParser;   
    DefaultNays=2;
    DefaultAngSpread=10;
%     DefaulRealization=1;
    DefaulRx=[1 1]; % user's array
    DefaulTxSector=[360 ,360];
    DefaulRxSector=[360 ,360];
    
    validScalarPosNum= @ (x) isnumeric(x) && isscalar(x) && (x > 0);
    ValidVectorPosNum=@(x) sum((isnumeric(x) & (~isscalar(x))& (x > 0))~=1)==0;
    validAngleSpread= @ (x) isnumeric(x) && isscalar(x) && (x >= 0);
    
    addRequired(p,'Tx',ValidVectorPosNum);
    addRequired(p,'K',validScalarPosNum);
    addRequired(p,'Ncl',validScalarPosNum);    
    addRequired(p,'F',validScalarPosNum);
    
    addParameter(p,'Rx',DefaulRx, ValidVectorPosNum);
    addParameter(p,'Nray',DefaultNays,validScalarPosNum);
    addParameter(p,'AngSpread',DefaultAngSpread,validAngleSpread);
    addParameter(p,'TxSector',DefaulTxSector,ValidVectorPosNum);
    addParameter(p,'RxSector',DefaulRxSector,ValidVectorPosNum);
%     addParameter(p,'Realization',DefaulRealization,validScalarPosNum);
%     addParameter(p,'Distribution',DefaulDistribution,@(x) any(validatestring(x,ExpectedDistribution)));
    parse(p,Tx,K,Ncl,F,varargin{:});
    
%     Realization=p.Results.Realization;
%     Distribution=p.Results.Distribution;
    Rx=p.Results.Rx;
    Nray=p.Results.Nray;
    AngSpread=p.Results.AngSpread *pi/180;
    TxSector=p.Results.TxSector*pi/180;
    RxSector=p.Results.RxSector*pi/180;
    F=p.Results.F;
%% Initial parameters    

Txz=Tx(1);
Txy=Tx(2);
if Txy==1
    error('BS can not be single antenna or An ULA on z-axis');
end

Nt=Txz*Txy;
BSArray='UPA';
if Txz==1
    BSArray='ULA';
end

Rxz=Rx(1);
Rxy=Rx(2);
if Rxz~=1 && Rxy==1
    error('User Antenna array must be on y-axis');
end
Nr=Rxz*Rxy;
UserArray='SingleAntenna';
if Rxy~=1
    UserArray='ULA';
    if Rxz~=1
        UserArray='UPA';
    end
end


H  = zeros(K,Nt,F);       % Channel matrix for realization repeat will be=  (Nr*Nt*realization)
At = zeros(Nt,Ncl*Nray,K); % At: transmit array response vectors, At is a Nt*(Nc*Nray) Matrix for ecah channel matrix
Ar = zeros(Nr,Ncl*Nray,K); % Ar: receive array response vectors, Ar is a Nr*(Nc*Nray) Matrix for ecah channel matrix
phi_t   = zeros(1,Ncl*Nray);         % The azimuth angle of AoDs
theta_t = zeros(1,Ncl*Nray);         % The elevation angle of AoDs
phi_r   = zeros(1,Ncl*Nray);         % The azimuth angle of AoAs
theta_r = zeros(1,Ncl*Nray);         % The elevation angle of AoAs
alpha   = zeros(Ncl*Nray,K);
L=Ncl*Nray;
gamma = sqrt((Nt*Nr)/(L));           %normalization factor

 %%   Path by Uniform distribution  
    if  Nray ~=1       
        for k = 1:K
            % Generate AoA & AOD for all rays
                % At first, we generate angle of each cluster with  the Uniform distribution
                AoD_Az = unifrnd(0,TxSector(1),Ncl,1); 
                AoD_El = unifrnd(0,TxSector(2),Ncl,1);
                AoA_Az = unifrnd(0,RxSector(1),Ncl,1);
                AoA_El = unifrnd(0,RxSector(2),Ncl,1);
            for c = 1:Ncl
                % at Second, we generate angle of rays with the Laplacian distribution  by using cluster angle as mean
                phi_t  ((c-1)*Nray+1 : Nray*c)   = Laprand(AoD_Az(c),AngSpread/sqrt(2),1,Nray);
                theta_t((c-1)*Nray+1 : Nray*c)   = Laprand(AoD_El(c),AngSpread/sqrt(2),1,Nray);
                phi_r  ((c-1)*Nray+1 : Nray*c)   = Laprand(AoA_Az(c),AngSpread/sqrt(2),1,Nray);
                theta_r((c-1)*Nray+1 : Nray*c)   = Laprand(AoA_El(c),AngSpread/sqrt(2),1,Nray);
            end 
            % If tranmitter array is ULA on y-axis
            if strcmp(BSArray,'ULA')
                theta_t(:)=pi/2;
            end
            % If user's array is ULA on y-axis
            if strcmp(UserArray,'ULA')
                theta_r(:)=pi/2;
            end
            % If users are single antenna
            if strcmp(UserArray,'SingleAntenna')
                theta_r(:)=pi/2;
                phi_r(:)=pi/2;
            end

            [Phi_t,~]=meshgrid(phi_t,(0:Nt-1));
            [Theta_t,~]=meshgrid(theta_t,(0:Nt-1));
            [Phi_r,~]=meshgrid(phi_r,(0:Nr-1));
            [Theta_r,~]=meshgrid(theta_r,(0:Nr-1));
            [~,TXy]=meshgrid(1:L,repmat((0:Txy-1),1,Txz));
            [~,TXz]=meshgrid(1:L,floor((0:Nt-1)/Txy));
            [~,RXy]=meshgrid(1:L,repmat((0:Rxy-1),1,Rxz));
            [~,RXz]=meshgrid(1:L,floor((0:Nr-1)/Rxy));
            At(:,:,k)=sqrt(1/Nt)*exp(1i* pi*(TXy.*sin(Phi_t).* sin(Theta_t)+ TXz.*cos(Theta_t) ));
            Ar(:,:,k)=sqrt(1/Nr)*exp(1i* pi*(RXy.*sin(Phi_r).* sin(Theta_r)+ RXz.*cos(Theta_r) ));
            alpha(:,k)= sqrt(1/2)*(randn(L,1) + 1j*randn(L,1));            
            Freq=(2*pi/F)*reshape(repmat([0:Ncl-1],Nray,1),[],1);
            for f = 1:F
                W=exp(-1j*Freq*f);
                Alpha(:,:,k,f)=gamma*diag(alpha(:,k).*W);
                H(k,:,f)=Ar(:,:,k)*Alpha(:,:,k,f)*At(:,:,k)';
            end        
        end
    end
% Reshape channel matrix
% if Nr==1
%     H=reshape(MuH,Nt,k).';
% else
%     H=MuH;
% end
end
