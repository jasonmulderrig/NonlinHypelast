function [P, R, K_T]=AssemblyNLElem(d,meshStruct)
% [P, R, K_T]=AssemblyNLElem(d,meshStruct);
% Calculate and assemble the global tangent stiffness matrix K_T, the 
% global internal force R vector, and the global external force P vector 
% using the global displacement vector d

% last update: 16 May 2021 C. Bonneville; J. Mulderrig; S. Srivatsa

% unpack necessary information
nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
numEls=meshStruct.numEls;
numEq=meshStruct.numEq;
gatherMat=meshStruct.gatherMat;

% initialize the global system
P = zeros(numEq,1);
R = zeros(numEq,1);
K_T = zeros(numEq,numEq);

% get the appropriate quadrature locations and weights
[qp,w] = Gauss(nnpe);
for e=1:numEls
    % Calculate local external force P_e vector, local internal force R_e
    % vector, and the local tangent stiffness matrix K_T_e for this element
    [P_e, R_e, K_T_e] = TwoDNLElem(d, meshStruct, e, qp, w);
    
    for Lrw = 1 : (nnpe*numDOF)
        Grw = gatherMat(e,Lrw); % global row index
        P(Grw) = P(Grw) + P_e(Lrw);
        R(Grw) = R(Grw) + R_e(Lrw);
        for Lcl = 1 : (nnpe*numDOF)
            Gcl = gatherMat(e,Lcl); % global column index
            
            % Assemble local tangent stiffness matrix into the global
            % tangent stiffness matrix here
            K_T(Grw,Gcl) = K_T(Grw,Gcl) + K_T_e(Lrw,Lcl);
        end
    end
end

% Package variables in global system as sparse objects to reduce solution
% time
P=sparse(P);
R=sparse(R);
K_T=sparse(K_T);
