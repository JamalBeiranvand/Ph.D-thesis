function [a,b] = findArea(width,varargin)
   defaultHeight = 1;
   defaultUnits = 'inches';
   defaultShape = 'rectangle';
   expectedShapes = {'square','rectangle','parallelogram'};

   p = inputParser;
   validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
   addRequired(p,'width',validScalarPosNum);
   addOptional(p,'height',defaultHeight,validScalarPosNum);
   addParameter(p,'units',defaultUnits,@isstring);
   addParameter(p,'shape',defaultShape,...
                 @(x) any(validatestring(x,expectedShapes)));
   parse(p,width,varargin{:});
   
   a = p.Results; 
   b=p;
end