function [Wopt,Fopt,Lambda]=DigitalBeamforming(H,Ns)
        [U,S,V] = svd(H);
        Fopt = V(:,1:Ns); 
        Wopt=U(:,1:Ns);
%           noisevar=1/SNR;
%         Wopt = ((Fopt'*(H'*H)*Fopt+noisevar*Ns*eye(Ns))\(Fopt'*H'))';
        Lambda=diag(S(1:Ns,1:Ns));
end