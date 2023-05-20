% This function is written according to
% "Doubling Phase Shifters for Efficient Hybrid Precoder Design in Millimeter-WaveCommunication Systems" 
% you can see the paper on this link: 
% 
% Inputs:
    % H:  Channel matrix
    % Ns=  Number of streams 
    % NRF: Nuber of RF chains
    % P: Power 
    % b: the number of bits in the resolution of phase shifters
% outputs
    % WRF: RF combiner matrix
    % WBB: baseband combiner matrix
    % FRF: RF precoder matrix
    % FBB: baseband precoder matrix
    
% Author  : Jamal Beiranvand (Jamalbeiranvand@gmail.com)
% google scholar: https://scholar.google.com/citations?user=S6LywwsAAAAJ&hl=en
%%
function [Wopt,Fopt]=FD_OFDM_MIMO(H,Ns,P)
    [Mr,Mt,F]=size(H);
    FOP=[];   
    for f=1:F
        [Wop,S,Fop]=svd(H(:,:,f));
        water_power = func_water_filling_algorithm((1)*ones(Ns,1),diag(S(1:Ns,1:Ns).^2),P);
        Fopt(:,:,f)=Fop(:,1:Ns)*sqrt(water_power);
        Wopt(:,:,f)=Wop(:,1:Ns);
    end
    
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function water_power = func_water_filling_algorithm(SNR,lambda,Trans_Power)

        Noise_Power= (1./(SNR.*lambda));
        Number_Channel= length(Noise_Power) ;
        [S_Number, dt]=sort(Noise_Power);
        sum(Noise_Power);
        for p=length(S_Number):-1:1
            T_P=(Trans_Power+sum(S_Number(1:p)))/p;
            Input_Power=T_P-S_Number;
            Pt=Input_Power(1:p);
            if(Pt(:)>=0)
                break
            end
        end
        Allocated_Power=zeros(1,Number_Channel);
        Allocated_Power(dt(1:p))=Pt;
        water_power = diag(Allocated_Power);
        end

end