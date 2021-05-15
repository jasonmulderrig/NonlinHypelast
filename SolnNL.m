function delta_d_plus_1 = SolnNL(G,K_T,boundStruct)
% delta_d_plus_1 = SolnNL(G,K_T,boundStruct);
% Solve the global matrix system of equations for the global displacement 
% increment vector delta_d_plus_1

% last update: 16 May 2021 C. Bonneville; J. Mulderrig; S. Srivatsa

% unpack the things you need
freeDOF = boundStruct.freeDOF;
numEq=length(G);

delta_d_plus_1 = zeros(numEq,1);
                             
K_T_FF	= K_T(freeDOF,freeDOF);    % Extract K_F matrix 
G_F = G(freeDOF);            % Extract G_F vector

delta_d_plus_1_F = -K_T_FF\G_F;

delta_d_plus_1(freeDOF) = delta_d_plus_1_F;





