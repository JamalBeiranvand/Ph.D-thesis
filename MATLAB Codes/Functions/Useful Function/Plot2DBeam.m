function Plot2DBeam(Input,ArrayStructTx,varargin)
p = inputParser;
defaultPhi=-180:180;
defaultTheta=0;
defaultArray='UPA';
% validScalarMatrix = @(x) isnumeric(x) && isscalar(x) && (x > 0);
expectedArrays = {'ULA','UPA'};
addParameter(p,'Array',defaultArray,@(x) any(validatestring(x,expectedArrays)));
addParameter(p,'Phi',defaultPhi);
addParameter(p,'Theta',defaultTheta);
parse(p,varargin{:});
Phi=p.Results.Phi;
Theta=p.Results.Theta;

if length(ArrayStructTx)==2
    Ntz=ArrayStructTx(1);
    Nty=ArrayStructTx(2);
    Nt=Ntz*Nty;
    Pt=floor((0:Nt-1)/Ntz);
    Qt=repmat((0:Ntz-1),1,Nty);
end
if length(ArrayStructTx)==1
    Nty=ArrayStructTx(1);
    Ntz=1;
    Nt=Ntz*Nty;
    Pt=floor((0:Nt-1)/Ntz);
    Qt=repmat((0:Ntz-1),1,Nty);
end
  if ( length(Phi)==1 || length(Theta)==1)
      fr=zeros(length(Phi),length(Theta));
      for i=1:length(Phi)
        for j=1:length(Theta)
            wave=sqrt(1/Nt)*exp(1i* pi*(Pt.*sind(Phi(i)).* cosd(Theta(j))-Qt.*sind(Theta(j)) ));        
            fr(i,j)=abs(sum(conj(wave)*Input));        
        end
      end        
  end

fr=fr/max(max(fr));
polarplot(Phi*pi/180,fr)
        
        
        