 
% At: transmit array response vectors, At is a Nt*(Nc*Nray) Matrix for ecah channel matrix
% Ar: receive array response vectors, Ar is a Nr*(Nc*Nray) Matrix for ecah channel matrix
% alpha: the gain of the lth ray in the ith propagation cluster
function [WRF,FRF]=AnalogBeamforming(At,Ar,alpha,Ns,Nq)
[~,num] = sort(abs(alpha),'descend');
FRF = At(:,num(1:Ns));
WRF = Ar(:,num(1:Ns));
      %% quantization
if (exist('Nq','var'))
    FRF=Qtz(FRF,Nq);
    WRF=Qtz(WRF,Nq);
end
