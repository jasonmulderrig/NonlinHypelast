function globalSystem=Assembly(MeshStruct)
% ASSEMBLY 
% Assemble global stiffness matrix K and global force vector F for the
% LinElast code. 

% last edit: 22 November 2017 H. Ritz

% unpack necessary information
nnpe=MeshStruct.nnpe;
numDOF=MeshStruct.numDOF;
numEls=MeshStruct.numEls;
numEq=MeshStruct.numEq;
gatherMat=MeshStruct.gatherMat;

% initialize the global system
Kt = zeros(numEq,numEq); 
F = zeros(numEq,1);
R = zeros(numEq,1);


% get the appropriate quadrature locations and weights
[qp,w] = Gauss(nnpe);
for e=1:numEls
    [kte,fe,re] = TwoDElemNL(MeshStruct, e, qp, w, d); % make the local stiffness matrix and 
                                              % local force vector for this element
    
    for Lrw = 1 : (nnpe*numDOF)
        Grw = gatherMat(e,Lrw); % global row index
        F(Grw) = F(Grw) + fe(Lrw);
        R(Grw) = R(Grw) + re(Lrw);
        for Lcl = 1 : (nnpe*numDOF)
            Gcl = gatherMat(e,Lcl); % global column index
            
            % Assemble local stiffness matrix into the global
            % stiffness matrix here
            Kt(Grw,Gcl) = Kt(Grw,Gcl) + kte(Lrw,Lcl);
        end
    end
end

% Package variables into the output structs. NOTE: both K and F are made
% sparse to greatly reduce solution time.
globalSystem.Kt=sparse(Kt);
globalSystem.F=sparse(F);
globalSystem.R=sparse(R);
