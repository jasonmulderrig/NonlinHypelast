% ------------------------------------------------------------------------|
%                                                                         |
% MAE4700-5700, Finite Element Analysis for Mechanical & Aerospace Design |
%                                                                         |
% Copyright: Cornell University (this software should not be used without |
% written permission)                                                     |
%                                                                         |
% Authors: N. Zabaras (zabaras@cornell.edu) & Xiang Ma (xm25@cornell.edu) |
% Last update: 17 Nov 2015 H. Ritz ; Y. Xu                                |
% ------------------------------------------------------------------------|
%
function  P=ApplyNaturalBC(i,boundStruct,meshStruct,P)
% Computed the element load from natural BCs.

% unpack necessary variables
sideInd = boundStruct.SurfNat(i,1);
qt = boundStruct.SurfNat(i,2); % tangetial traction
qn = boundStruct.SurfNat(i,3); % normal traction
BElems = boundStruct.elements(sideInd).Elems;
SurfID = boundStruct.elements(sideInd).SurfaceIndicator;

nnpe=meshStruct.nnpe;
numDOF=meshStruct.numDOF;
nCoords=meshStruct.nCoords;
elCon=meshStruct.elCon;
gatherMat=meshStruct.gatherMat;
ldof = nnpe*numDOF;       % Degrees of freedom per element

Num = length(BElems); %  Extract the number of boundary elements

pt=1/sqrt(3);          % Set up Gauss quadrature for boundary integral
gp = [-pt,pt];         % 2 points for one dimension - this will need to
w  = [1,1];            % be generalized for other situations.

for elmID = 1 : Num
    
    P_e = zeros(ldof,1);  % initialize elemental external work vector
    
    glb = elCon(BElems(elmID),:);   % Global number of the element nodes
    
    coord = nCoords(glb,:);          % Global nodal coordinates for the element
    % nodes as a column vector
    
    for qp = 1:length(w)             % loop over all the gauss points on the
        % boundary segment
        
        %--------------------------------------------------------------------
        % First set up the information of the finite element for the current
        % integration point. In other words, put the value and derivative of
        % basis function, the physical coordinates of the integration points in
        % the struct fe.
        %------------------------------------------------------------------
        %
        % Different basis functions according to different sides
        
        if ( nnpe == 4)
            switch (SurfID(elmID))
                case -2
                    n = [(1 - gp(qp))/2, (1 + gp(qp))/2, 0, 0];
                    dn = [-1/2, 1/2, 0, 0];
                case +1
                    n = [0, (1 - gp(qp))/2, (1 + gp(qp))/2, 0];
                    dn = [0, -1/2, 1/2, 0];
                case +2
                    n = [0, 0, (1 - gp(qp))/2, (1 + gp(qp))/2];
                    dn = [0, 0, -1/2, 1/2];
                case -1
                    n = [(1 + gp(qp))/2, 0, 0, (1 - gp(qp))/2];
                    dn = [1/2, 0, 0, -1/2];
            end
            
        elseif  ( nnpe == 3)
            switch (SurfID(elmID))
                case -2
                    n = [(1 - gp(qp))/2, (1 + gp(qp))/2, 0];
                    dn = [-1/2, 1/2, 0];
                case +1
                    n = [0, (1 - gp(qp))/2, (1 + gp(qp))/2];
                    dn = [0, -1/2, 1/2];
                case -1
                    n = [(1 + gp(qp))/2, 0, (1 - gp(qp))/2];
                    dn = [1/2, 0, -1/2];
            end
        else
            fprintf(1, 'This element type is not implemented\n');
        end
        
        P_elm.N = n;           % put basis function into struct P_elem
        P_elm.x = n*coord;     % Calculate the global coordinates of the
        % current integration points
        Jacc = dn*coord;      % Calculate the components of the Jaccobian matrix
        % Now it is a vector
        detJ = norm(Jacc);    % Calculate the norm of the vector
        
        P_elm.detJxW = detJ * w(qp); % Calculate the integration weight
        
        P_elm.nx = Jacc(2)/detJ;    % Calculate the direction cosines
        P_elm.ny = -Jacc(1)/detJ;

        % Extract the N matrix
        if (nnpe == 4)
            N = [ P_elm.N(1),   0,    P_elm.N(2),   0,     P_elm.N(3),     0,  P_elm.N(4),  0;
                0 ,  P_elm.N(1),    0   ,  P_elm.N(2),   0,    P_elm.N(3),  0 ,  P_elm.N(4)];
        else
            N = [ P_elm.N(1),   0,    P_elm.N(2),   0,     P_elm.N(3),     0;
                0 ,  P_elm.N(1),    0   ,  P_elm.N(2),   0,    P_elm.N(3)];
            
        end
        
        % Transform to x and y direction using the following formula
        % where nx and ny are the direction cosines at the current boundary
        q = [P_elm.nx -P_elm.ny;P_elm.ny  P_elm.nx] * [qn;qt];
        
        % Assemble into force vector
        P_e = P_e +  N' * q * P_elm.detJxW;
        
        
    end
    clear P_elm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble to global
    glbROW = gatherMat(BElems(elmID),:); % global row index
    P(glbROW)=P(glbROW)+P_e;

end