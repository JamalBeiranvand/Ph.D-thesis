function a=Testfunc(x,varargin)
expectedShapes = {'square','rectangle','parallelogram','triangle'};
p = inputParser;
defaultShape='square';
addParameter(p,'shape',defaultShape,@(x) any(validatestring(x,expectedShapes)));
% addRequired(p,'x');
parse(p,varargin{:});
shapeName=p.Results.shape;
switch shapeName
        case {'square','rectangle'}
            a = x;
        case {'triangle'}
            a = x^2;
        otherwise
            error('Unknown shape passing validation.')
end