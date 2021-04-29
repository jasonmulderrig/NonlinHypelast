% GlobalSystem()
% This function determines the global residual vector G and global 
% tangential stiffness vector KT for each iteration, and applies both
% essential and natural BCs to the global system
% NOTE: FUNCTION CALLS DO NOT INCLUDE INPUTS/OUTPUTS IN FINAL FORM
function [G,KT] = GlobalSystem(vi,loadstep)

% Calculate global residual vector G, global tangential stiffness vector KT
[G,KT] = Assembly();

% Calculate and apply natural, essential BCs to global system of equations