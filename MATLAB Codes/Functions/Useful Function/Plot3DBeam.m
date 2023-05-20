function [X,Y,Z]=Plot3DBeam(Input,ArrayStruct,varargin)
p = inputParser;
defaultPhi=-180:180;
defaultTheta=-90:90;
defaultArray='UPA';
% validScalarMatrix = @(x) isnumeric(x) && isscalar(x) && (x > 0);
expectedArrays = {'ULA','UPA'};
addParameter(p,'Array',defaultArray,@(x) any(validatestring(x,expectedArrays)));
addParameter(p,'Phi',defaultPhi);
addParameter(p,'Theta',defaultTheta);
parse(p,varargin{:});
Phi=p.Results.Phi;
Theta=p.Results.Theta;

if length(ArrayStruct)==2
    Ntz=ArrayStruct(1);
    Nty=ArrayStruct(2);
    Nt=Ntz*Nty;
    Pt=floor((0:Nt-1)/Ntz);
    Qt=repmat((0:Ntz-1),1,Nty);
end
if length(ArrayStruct)==1
    Nty=ArrayStruct(1);
    Ntz=1;
    Nt=Ntz*Nty;
    Pt=floor((0:Nt-1)/Ntz);
    Qt=repmat((0:Ntz-1),1,Nty);
end

 if ( length(Phi)>1 && length(Theta)>1)
        [THETA,PHI]=meshgrid(Theta,Phi);
        fr=zeros(size(THETA));
        for i=1:length(Phi)
            for j=1:length(Theta)
                wave=sqrt(1/Nt)*exp(1i* pi*(Pt.*sind(Phi(i)).* cosd(Theta(j))-Qt.*sind(Theta(j)) ));        
                fr(i,j)=abs(sum(conj(wave)*Input));        
            end
        end
        fr=fr/max(max(fr));

        %figure;

        hold on
        [X,Y,Z] = sph2cart(PHI*pi/180,THETA*pi/180,fr);
        C=abs(fr);
        surf(X,Y,Z,C);
        ax=gca;
        ax.XLim=[min(min(X)) max(max(X))];
        ax.YLim=[min(min(Y)) max(max(Y))];
        ax.ZLim=[min(min(Z)) max(max(Z))];
        hold on
% 
%         C=abs(X.^2+Y.^2);
%         z=Z(:,:);
%         z(:,:)=ax.ZLim(1);
%         surf(X,Y,z,C)
% 
%         y=Y(:,:);
%         y(:,:)=ax.YLim(2);
%         C=abs(Z.^2+X.^2);
%         surf(X,y,Z,C)


        shading interp 
        lighting gouraud 
        material dull
        camlight('headlight')

%         xlabel('x');
%         ylabel('Y');
%         grid on
%         box on
%         hold on
%         plot3([ax.XLim(1) ax.XLim(2)],[0 0],[ax.ZLim(1) ax.ZLim(1)], 'Color', [.6 .2 .2])
%         plot3([0 0],[ax.YLim(1) ax.YLim(2)],[ax.ZLim(1) ax.ZLim(1)], 'Color', [.6 .2 .2])
%         plot3([0 0],[ax.YLim(2) ax.YLim(2)],[ax.ZLim(1) ax.ZLim(2)],'Color', [.6 .2 .2])
%         plot3([ax.XLim(1) ax.XLim(2)],[ax.YLim(2) ax.YLim(2)],[0 0],'Color', [.6 .2 .2])
%         ax.XTick = [-1 0 1];
%         ax.YTick = [-1 0 1];
%         ax.ZTick = [-1 0 1];
%         ax.XLabel.String='x';
%         ax.YLabel.String='y';
%         ax.ZLabel.String='z';
        set(gcf,'color','w')
        view([ -1 ,-1,0]);
        ax.TickLabelInterpreter='latex';
        ax.YLabel.Interpreter='latex';
        ax.XLabel.Interpreter='latex';
        ax.FontSize = 11;
        ax.Visible='off';
 end