% GlobalSystemCalcn()
% This function determines the global residual vector G and global 
% tangential stiffness vector KT for each iteration, and applies both
% essential and natural BCs to the global system
% NOTE: FUNCTION CALLS DO NOT INCLUDE INPUTS/OUTPUTS IN FINAL FORM
function [G,K_T] = GlobalSystemCalcn(d,meshStruct,boundStruct)

% Calculate global residual vector G, global tangential stiffness vector KT
[P, R, K_T] = AssemblyNLElem(d,meshStruct);

% Calculate and apply natural BCs to global system of equations

% Apply natural BC(s)
for i = 1:size(boundStruct.SurfNat,1)    % Loop over all natural boundaries
    P = ApplyNaturalBC(i,boundStruct,meshStruct,P); 
end

G = R-P;