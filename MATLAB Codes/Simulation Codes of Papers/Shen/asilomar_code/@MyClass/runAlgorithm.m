function [ schedule, V ] = runAlgorithm( obj, weight )

numIter = 50;

switch obj.algorithm
    case 'proposed'
        [ schedule, V ] = runProposed(obj, weight, numIter);
    case 'WMMSE'
        [ schedule, ~, V ] = runWMMSE(obj, weight, numIter);
    case 'fixed'
        [ schedule, ~, V ] = runFixed(obj, weight, numIter);
    case 'Beam'   
        [ schedule, ~, V ] = runBeam( obj, weight );
    otherwise
        error('unidentified algorithm')
end

end