clear
clc
 addpath(genpath('C:\Users\Jamal\OneDrive\Jamal\MATLAB Code\Functions' ))% Add path of the "Functions" folder.
N=17;
[Setb,Sb,Set,S]=GenUniqueSet(N);
%scatterplot(Set)
Eta=max(abs(Set));

NData=1e4;
% r=rand(NData,1);
% t=2*pi*rand(NData,1);
% z=r.*exp(1j*t);
% scatterplot(z)
% 
%z=randn(NData,1)+1j*randn(NData,1);
% scatterplot(z)

z=-Eta+2*Eta*rand(NData,1)+1j*(-Eta+2*Eta*rand(NData,1));
z=z(abs(z)<Eta);
%scatterplot(z)

Zm=MUMISO_DownModified(z,Set,S);
% hold on
% scatterplot(Zm)

MSe=mean(abs(Zm-z).^2)/Eta

