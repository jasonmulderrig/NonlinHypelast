% GlobalSystemCalcn()
% This function determines the global residual vector G and global 
% tangential stiffness vector KT for each iteration, and applies both
% essential and natural BCs to the global system
% NOTE: FUNCTION CALLS DO NOT INCLUDE INPUTS/OUTPUTS IN FINAL FORM
function [G,KT] = GlobalSystemCalcn(d,meshStruct,boundStruct)

% Calculate global residual vector G, global tangential stiffness vector KT
[P, R, K_T] = AssemblyNLElem(d,meshStruct);

% Calculate and apply natural, essential BCs to global system of equations