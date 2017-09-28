function coOrds3d = CoordinateTo3d(coOrds2d,isLeft)
%-------------------------------------------------------------------------------
% Uses formulation of Varier et al. to infer a 3d coordinate from the planar 2d
% coordinate: C.Pos;
% [C.RegionStruct.isLeft]
% EXAMPLE: coOrds3d = CoordinateTo3d(C.Pos,[C.RegionStruct.isLeft]);
%-------------------------------------------------------------------------------

radiusWorm = 50e-3; % (mm); 50 micron = 50e-3 mm.

numNeurons = size(coOrds2d,1);
coOrds3d = zeros(numNeurons,3);
coOrds3d(:,1:2) = coOrds2d;
% Put neurons on the circumpherence of the cross section using Pythagorus:
coOrds3d(:,3) = (sqrt(radiusWorm^2 - coOrds2d(:,2).^2))/2;
% Swap left/right in z for bilerateral neurons:
coOrds3d(isLeft,3) = - coOrds3d(isLeft,3);

end
