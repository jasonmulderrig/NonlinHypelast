function [G,K_T] = GlobalSystemCalcn(d,meshStruct,boundStruct)
% [G,K_T] = GlobalSystemCalcn(d,meshStruct,boundStruct);
% Calculate the global residual G vector and the global tangent stiffness 
% matrix K_T using the displacement vector d

% last update: 16 May 2021 C. Bonneville; J. Mulderrig; S. Srivatsa

% Calculate global residual vector G, global tangential stiffness vector KT
[P, R, K_T] = AssemblyNLElem(d,meshStruct);

% Apply natural BCs to global system of equations
for i = 1:size(boundStruct.SurfNat,1) % Loop over all natural boundaries
    P = ApplyNaturalBC(i,boundStruct,meshStruct,P); 
end

G = R-P;