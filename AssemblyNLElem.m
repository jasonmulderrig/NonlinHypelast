function [P, R, K_T]=AssemblyNLElem(d,meshStruct)
% ASSEMBLY 
% Assemble global stiffness matrix K and global force vector F for the
% LinElast code. 

% last edit: 22 November 2017 H. Ritz

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
    [P_e, R_e, K_T_e] = TwoDNLElem(d, meshStruct, e, qp, w); % make the local stiffness matrix and 
                                              % local force vector for this element
    
    for Lrw = 1 : (nnpe*numDOF)
        Grw = gatherMat(e,Lrw); % global row index
        P(Grw) = P(Grw) + P_e(Lrw);
        R(Grw) = R(Grw) + R_e(Lrw);
        for Lcl = 1 : (nnpe*numDOF)
            Gcl = gatherMat(e,Lcl); % global column index
            
            % Assemble local stiffness matrix into the global
            % stiffness matrix here
            K_T(Grw,Gcl) = K_T(Grw,Gcl) + K_T_e(Lrw,Lcl);
        end
    end
end

% Package variables into the output structs. NOTE: both K and F are made
% sparse to greatly reduce solution time.
P=sparse(P);
R=sparse(R);
K_T=sparse(K_T);
