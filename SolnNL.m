% SOLN 
% Apply the essential boundary conditions and solve the global system for
% the nodal degrees of freedom for LinElast. 
% last edit: 22 November 2017 H. Ritz

function delta_d_plus_1 = SolnNL(G,K_T,boundStruct)
% unpack the things you need
essDOF = boundStruct.essDOF;
numEq=length(G);

freeDOF=setdiff(1:numEq,essDOF); % this returns the indices to the DOF that 
                             % do NOT have essential boundary conditions 
                             % (free DOF)

delta_d_plus_1 = zeros(numEq,1);
                             
K_T_FF	= K_T(freeDOF,freeDOF);    % Extract K_F matrix 
G_F = G(freeDOF);            % Extract G_F vector

delta_d_plus_1_F = K_T_FF\( -G_F );
 
delta_d_plus_1(freeDOF) = delta_d_plus_1_F;





