%  clear
%  clc
%  warning off
function [F,VT,Performance_opt]=Mainfunction_UpdateF(HT,NRF,P)

%Initialization
%Set parameters Nt, K, Nrf, and P    

    Q=size(HT,2);
    Nt=size(HT,2);
    K=size(HT,3);
    N=size(HT,4);
    T=NRF;
    for f=1:N
        for k=1:K
            Ht(k,:,f)=HT(1,:,k,f);
        end
    end

    Iterations= 20; % the number of iteration
    R=0; % R is a local variable to compute performance
    Y=zeros(K,N);
    Gamma=zeros(K,N);
    [~,~,Set,~]=GenUniqueSet(11);
%Generate H[n] for all n; Ht includes all H[n];
  %  Ht=sqrt(1/2)*(randn(K,Nt,N)+1j*randn(K,Nt,N));
%Calculate the optimal performance from  the function Chen_Alg5
    % the optimal performance
    [Performance_opt,Fopt]=Chen_Alg5(Ht,P);
    disp('The optimal Performance is =')
    disp(Performance_opt)
%Initialize F and V for all n
   [F,Vt]=InitializeFV(Ht,Fopt,NRF,P); %Vt includes all V[n]
%% Iterative Algorithm 
%  for iter = 1:Iterations % 
sum_rate = 10e2;
epsilon =N*10^(-3);
R = 0;
 while abs(sum_rate - R) >= epsilon 
     sum_rate = R;
         % Step 1: 
       % Fixed F, V, and Y, Update gamma for all n 
       for n=1:N 
           H=Ht(:,:,n); % select H number n from Ht  
           V=Vt(:,:,n); % select V number n from Vt 
           Gamma(:,n) = updateGamma(H,F,V);
       end
       
        %Step 2: 
       %% Fixed F, V, and gamma , Update y for all n
       for n=1:N 
           H=Ht(:,:,n); % select H number n from Ht 
           V=Vt(:,:,n); % select V number n from Vt 
           Y(:,n) = updateY(H,F,V, Gamma(:,n));
       end
      
        % Step 3: 
        %% update F by fixing other parameters
        % here we have to add a function to update F
%         tic
        F=updateF(F,Ht,K,N,Nt,NRF,Gamma,Vt,Y,P);
%        toc
        % Step 4: 
        %% Fixed F, gamma, and Y, Update V for all n 
        for n=1:N
            H=Ht(:,:,n); % select H number n from Ht 
            Vt(:,:,n)  = updateV(H, F,Y(:,n), Gamma(:,n),P);        
        end   
        
       % check performance
        R = 0;
        for n=1:N
            H=Ht(:,:,n); % select H number n from Ht 
            V=Vt(:,:,n); % select V number n from Vt 
            R=R+Rates2(H,F*V);
        end 
       % disp(R/N) ;
        
        % Non-decreasing
        if R < sum_rate
            R = sum_rate;
            break; % end loop while-do
        end
 end
 for f=1:N
        for k=1:K
            VT(:,1,k,f)=Vt(:,k,f);
        end        
 end
 
 disp("=== End of Problem ===");
end
 

%% Functions

%Update F

%Version 1:
function F = updateF(F,Ht,K,N,Nt,NRF,Gamma,Vt,Y,P)
    % So I expect to get a column vector g: 1st NRF elements - Real part, next NRF - Imag
    % Initial vector g
    a = F;
    a1 = [];
    for i = 1:Nt
        for j = 1:NRF        
            a1 = [a1 a(i,j)]; % a1 = a1_r + 1i*a1_i;   
        end
    end
    a1_r = real(a1);
    a1_i = imag(a1);
    
    a11 = [a1_r a1_i];
    
    % Main function
    objective = @(g) f_objective(g, F, Ht, Vt, Nt, NRF, K, N, Gamma);
    
    % x0 = ones(1,2*Nt*NRF);
    x0 = a11;
    options = optimoptions('fmincon','MaxFunctionEvaluations',3000);
    [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]...
        = fmincon(objective,x0,[],[],[],[],[],[],@(g)nonlincon(g, N, Nt, NRF, Vt, P),options);
    
    % Convert F matrix
    F_out = [];
    a2_r = [];
    for ii = 1:NRF:Nt*NRF
        a2_r = [a2_r; X(1,[ii:(ii+NRF-1)])];
    end
    
    a2_i = [];
    for iii = (Nt*NRF+1):NRF:(2*Nt*NRF)
        a2_i = [a2_i; X(1,[iii:(iii+NRF-1)])];
    end
    
    F_out = a2_r + 1i*a2_i; % convert to complex matric
    
    F = F_out;
end

function [object] = f_objective(g, F, Ht, Vt, Nt, NRF, K, N, Gamma)
    a2 = [];
    a2_r = [];
    for ii = 1:NRF:Nt*NRF
        a2_r = [a2_r; g(1,[ii:(ii+NRF-1)])];
    end
    
    a2_i = [];
    for iii = (Nt*NRF+1):NRF:(2*Nt*NRF)
        a2_i = [a2_i; g(1,[iii:(iii+NRF-1)])];
    end
    
    a2 = a2_r + 1i*a2_i; % convert to complex matric F
    
    % Given other variables
    delta11 = nan(K,N);
    delta21 = nan(K,N);
    delta_sum  = nan(K,N);
    for n=1:N
        H=Ht(:,:,n); % select H number n from Ht
        V=Vt(:,:,n); % select V number n from Vt
        %            Y(:,n) = updateY(H,F,V, Gamma(:,n));
        
        V1 = F*V; %consider V=F*Vn for nth frequence of previous iteration
        [K,~]=size(H);
        y = nan(K,1);
        delta1 = nan(K,1);
        delta2 = nan(K,1);
        for i = 1:K
            hi=H(i,:);
            Vi=V1(:,i);
            Vi_opt = V(:,i);
            B = 1;
            B2 = 1;
            for j = 1:K
                Vj=V1(:,j);
                Vj_opt = V(:,j);
                B = B + (hi*Vj*(Vj')*hi');
                B2 = B2 + (hi*a2*Vj_opt*(Vj_opt')*(a2)'*hi');
            end
            % Calculate y based on the F of previous iteration
            y(i) = (sqrt(1+Gamma(i,n)) *hi*Vi)/B;
            % delta1
            delta1(i) = sqrt(1+Gamma(i,n))*2*real(y(i).'*hi*a2*Vi_opt);
            delta2(i) = y(i).'*B2*y(i);
        end
    %     Y(:,n) = y;
        delta11(:,n) = delta1;
        delta21(:,n) = delta2;
        delta_sum(:,n) = delta1 - delta2;
    end
    % Output
    object = - real(sum(delta_sum(:)));
end

function [c,ceq] = nonlincon(g, N, Nt, NRF, Vt, P)
    a2 = [];
    a2_r = [];
    for ii = 1:NRF:Nt*NRF
        a2_r = [a2_r; g(1,[ii:(ii+NRF-1)])];
    end
    
    a2_i = [];
    for iii = (Nt*NRF+1):NRF:(2*Nt*NRF)
        a2_i = [a2_i; g(1,[iii:(iii+NRF-1)])];
    end
    
    a2 = a2_r + 1i*a2_i; % convert to complex matric
    
    for n=1:N
        V=Vt(:,:,n); 
        c(n) = real(trace(a2*V*(V')*a2')) - P - 10^(-4);
    end
ceq = [];
end

%% Initialization Step
function [F,V]=InitializeFV(Ht,Fopt,T,P)
 [K,~,N]=size(Ht);
  
    FOP=[];
    for k=1:K
        for f=1:N
            FOP=[FOP,Fopt(:,k,f)];
        end
    end
    Ns=1;
    [UF,S,VF]=svd(FOP); 
    FRF=UF(:,1:T);
    Fbb=S(1:T,1:T)*VF(:,1:T)';
    for k=1:K
        for f=1:N
            SColumn=(k-1)*N*Ns+(f-1)*Ns+1;
            EColumn=(k-1)*N*Ns+f*Ns;
            FBB(:,k,f)=Fbb(:,SColumn:EColumn);  
            FBB(:,k,f)=sqrt(Ns*P)*FBB(:,k,f)/norm(FRF*FBB(:,k,f),'fro');
        end
    end
    
    F=FRF;
    V=squeeze(FBB);
end
%% Gamma
function [ gamma ] = updateGamma(H,F,Vn) 
V=F*Vn; %consider V=F*Vn for nth frequence
[K,~]=size(H);
gamma = zeros(K,1);

    for i = 1:K
        hi=H(i,:);
        Vi=V(:,i);
        B = 1;
        for j= 1:K
            Vj=V(:,j);
            B = B + hi*Vj*(Vj)'*hi';
        end  
        B=B-hi*Vi*(Vi)'*hi';
        gamma(i) = (Vi)'*hi'*(B\(hi*Vi));
    end
end
%% Y
function [ y ] = updateY(H, F, Vn, gamma)
V=F*Vn; %consider V=F*Vn for nth frequence
[K,~]=size(H);
    y = nan(K,1);
    for i = 1:K
        hi=H(i,:);
        Vi=V(:,i);
        B = 1;
        for j = 1:K
            Vj=V(:,j);
            B = B + (hi*Vj*(Vj')*hi');
        end  
        y(i) = B\(sqrt(1+gamma(i)) *hi*Vi);
%         if isnan(y(i))
%             y(i) = 0;
%         end
    end
end
%% V
function [ V ] = updateV(Hn, F,y, r,P)
H=Hn*F;
[K,NRF]=size(H);
    % update V       
        G= zeros(NRF,NRF);
        for k = 1:K
            G = G+ H(k,:)'*y(k)*y(k)'*H(k,:);
        end 
        %G=F'*G*F;
        
        tempV=zeros(NRF,K);
        for k=1:K
            A = sqrt(1+r(k))*H(k,:)'*y(k);
            tempV(:,k) = G\A;
        end
        if trace(F*tempV*(tempV')*F')<=P 
            V=tempV;
            V=sqrt(P)*V/sqrt(trace(F*V*(V')*F'));
        else
            %% Compute the Optimal mu
            Mustar=Bisection(H,F,y,r,P);
            T=G+Mustar*(F')*F;
            for k=1:K
                A = sqrt(1+r(k))*H(k,:)'*y(k);
                V(:,k) = T\A;
            end 
           
        end  
end

%% Using Bisection method to find the optimal mu
function Mu_opt=Bisection(H,F,y,r,P)
       [K,NRF]=size(H);
       B = zeros(NRF,NRF);
        for k = 1:K
            B = B+ H(k,:)'*y(k)*y(k)'*H(k,:);
        end 
        LRight = 1;
        while 1
            J=B+LRight*(F')*F;
            for k=1:K
                A = sqrt(1+r(k))*H(k,:)'*y(k);
                tempV(:,k)= J\A;
            end
            Ptxf=trace(F*tempV*(tempV')*F');
            if Ptxf < P
                break                     
            end  
            LRight = 5*LRight;
        end
        LLeft=0;
        while 1
            Lambda=mean([LLeft,LRight]);
            J=B+Lambda*(F')*F;
            for k=1:K
                A = sqrt(1+r(k))*H(k,:)'*y(k);
                tempV(:,k)= J\A;
            end
            Ptxf=trace(F*tempV*(tempV')*F');
            if abs(Ptxf-P) < P/1e4
                Mu_opt= Lambda;
                break 
            end
            if Ptxf > P
                LLeft = Lambda;
            else
                LRight = Lambda;
            end
            
            if abs(LRight-LLeft)<10e-10
                Mu_opt= Lambda;
                break
            end
        end         
            
end

%% Sum Rate
function R=Rates2(H,FBB)
[KK,~]=size(H);
Rfk=zeros(1,KK);
for k=1:KK
  Fkf=FBB(:,k);
  Hkf= H(k,:);
  
  OMEGA=0;
  for kp=1:KK
      Fkp=FBB(:,kp);
      OMEGA=OMEGA+Hkf*Fkp*(Fkp')*Hkf';
  end               
  Omega=OMEGA-Hkf*Fkf*(Fkf')*Hkf';
  Omega=(Omega+1);              
  Rfk(k)=log2(det(1 + (Fkf)'  * (Hkf') *(Omega\(Hkf * Fkf))));
end          
     R=sum(Rfk);
end


