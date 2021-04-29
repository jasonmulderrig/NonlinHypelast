function [kte, fe, re] = TwoDElemNL(meshStruct, elmID, qp, w, d)
% [localstiffnessmatrix, localforcevector] = TwoDElem(meshStruct,elementnumber,QuadPoints,Weights)
% generate the local stiffness matrix and local force vector from
% body forces for use with LINELAST code.
% last edit: 1 May 2015 H. Ritz

% unpack necessary information
nnpe = meshStruct.nnpe;
numDOF = meshStruct.numDOF;
nCoords = meshStruct.nCoords;
elCon = meshStruct.elCon;
mu = meshStruct.Material.mu;
lambda = meshStruct.Material.lambda;

D = [lambda + 2 * mu, lambda, 0;
     lambda, lambda + 2 * mu, 0;
     0, 0, mu]; % St Venant formulation

kte = zeros(numDOF * nnpe);   % initialize elemental stiffness matrix
fe = zeros(numDOF * nnpe, 1); % initialize elemental force vector
re = zeros(numDOF * nnpe, 1); % initialize elemental residual vector

% get nodal coordinates for this element
xvec = nCoords(elCon(elmID,:), 1); % these are column vectors of the
yvec = nCoords(elCon(elmID,:), 2); % physical coordinates of the
                                % nodes of this element

ux = xvec + d(elCon(elmID, :) * 2 - 1);                  
uy = yvec + d(elCon(elmID, :) * 2);

for iqp = 1 : size(qp, 1) % loop over quadrature points
    % for each quadrature point ...
    NofXiEta = Nmatrix(qp(iqp, :), nnpe); % row vector of shape functions
                                      % at this QP in parent domain
    dNdXiEta = dNmatrix(qp(iqp, :),nnpe);% 2xnnpe array of shape func derivs
                                      % at this QP in parent domain
    JofXiEta = dNdXiEta * [xvec, yvec];    % jacobian at this QP (2x2 matrix)
    dNdXY = inv(JofXiEta) * dNdXiEta;     % 2xnnpe array of dNdX at ths QP
    
    
    % --------- Deformation gradient Fe (Slide 13/18 - Eq 24) ----------- %
    
    F1 = [ux(1) * dNdXY(1, 1), ux(1) * dNdXY(2, 1);
          uy(1) * dNdXY(1, 1), uy(1) * dNdXY(2, 1)]; % Fe for node 1
      
    F2 = [ux(2) * dNdXY(1, 2), ux(2) * dNdXY(2, 2);
          uy(2) * dNdXY(1, 2), uy(2) * dNdXY(2, 2)]; % Fe for node 2
      
    F3 = [ux(3) * dNdXY(1, 3), ux(3) * dNdXY(2, 3);
          uy(3) * dNdXY(1, 3), uy(3) * dNdXY(2, 3)]; % Fe for node 3
      
    F4 = [ux(4) * dNdXY(1, 4), ux(4) * dNdXY(2, 4);
          uy(4) * dNdXY(1, 4), uy(4) * dNdXY(2, 4)]; % Fe for node 4
      
    Fe = F1 + F2 + F3 + F4; % Deformation gradient for given QP (Sum or Fe_k)
    
    % --------- Euler Lagrange Tensor Ee (Slide 13/18 - Eq 25) ---------- %
    
    Ee = 0.5 * (Fe' * Fe - eye(2));
    
    % ------------------ BL Matrix (Slide 14/18 - Eq 26) ---------------- %
    
    BL1 = [Fe(1, 1) * dNdXY(1, 1), Fe(2, 1) * dNdXY(1, 1);
           Fe(1, 2) * dNdXY(2, 1), Fe(2, 2) * dNdXY(2, 1);
           Fe(1, 1) * dNdXY(2, 1) + Fe(1, 2) * dNdXY(1, 1),  Fe(2, 1) * dNdXY(2, 1) + Fe(2, 2) * dNdXY(1, 1)]; % BL at node 1
    
    BL2 = [Fe(1, 1) * dNdXY(1, 2), Fe(2, 1) * dNdXY(1, 2);
           Fe(1, 2) * dNdXY(2, 2), Fe(2, 2) * dNdXY(2, 2);
           Fe(1, 1) * dNdXY(2, 2) + Fe(1, 2) * dNdXY(1, 2),  Fe(2, 1) * dNdXY(2, 2) + Fe(2, 2) * dNdXY(1, 2)]; % BL at node 2
       
    BL3 = [Fe(1, 1) * dNdXY(1, 3), Fe(2, 1) * dNdXY(1, 3);
           Fe(1, 2) * dNdXY(2, 3), Fe(2, 2) * dNdXY(2, 3);
           Fe(1, 1) * dNdXY(2, 3) + Fe(1, 2) * dNdXY(1, 3),  Fe(2, 1) * dNdXY(2, 3) + Fe(2, 2) * dNdXY(1, 3)]; % BL at node 3
    
    BL4 = [Fe(1, 1) * dNdXY(1, 4), Fe(2, 1) * dNdXY(1, 4);
           Fe(1, 2) * dNdXY(2, 4), Fe(2, 2) * dNdXY(2, 4);
           Fe(1, 1) * dNdXY(2, 4) + Fe(1, 2) * dNdXY(1, 4),  Fe(2, 1) * dNdXY(2, 4) + Fe(2, 2) * dNdXY(1, 4)]; % BL at node 4
       
    % ------------------ PK Stress (Slide 12/18 - Eq 23) ---------------- %
    
    Se = D * [Ee(1, 1); Ee(2, 2); 2 * Ee(1, 2)]; % Piola-Kirchoff stress tensor in voight notation for the element
    
    % ------------------ G Matrix (Slide 16/18 - Eq 31) ---------------- %
    
    Se_matrix = [Se(1), Se(3); Se(3), Se(2)]; % Piola-Kirchoff stress tensor in matrix notation for the element
    
    G11 = dNdXY(:, 1)' * Se_matrix * dNdXY(:, 1);
    G12 = dNdXY(:, 1)' * Se_matrix * dNdXY(:, 2);
    G13 = dNdXY(:, 1)' * Se_matrix * dNdXY(:, 3);
    G14 = dNdXY(:, 1)' * Se_matrix * dNdXY(:, 4);
    
    G21 = dNdXY(:, 2)' * Se_matrix * dNdXY(:, 1);
    G22 = dNdXY(:, 2)' * Se_matrix * dNdXY(:, 2);
    G23 = dNdXY(:, 2)' * Se_matrix * dNdXY(:, 3);
    G24 = dNdXY(:, 2)' * Se_matrix * dNdXY(:, 4);
    
    G31 = dNdXY(:, 3)' * Se_matrix * dNdXY(:, 1);
    G32 = dNdXY(:, 3)' * Se_matrix * dNdXY(:, 2);
    G33 = dNdXY(:, 3)' * Se_matrix * dNdXY(:, 3);
    G34 = dNdXY(:, 3)' * Se_matrix * dNdXY(:, 4);
    
    G41 = dNdXY(:, 4)' * Se_matrix * dNdXY(:, 1);
    G42 = dNdXY(:, 4)' * Se_matrix * dNdXY(:, 2);
    G43 = dNdXY(:, 4)' * Se_matrix * dNdXY(:, 3);
    G44 = dNdXY(:, 4)' * Se_matrix * dNdXY(:, 4);
    
    XY = NofXiEta * [xvec, yvec]; % physical coordinate of this QP (1x2)
    
    % now create the N matrix
    N = []
    for np=1:nnpe
        N = [N,[NofXiEta(np), 0; 0, NofXiEta(np)]];
    end
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental body force vector
    
    kte_qp = [BL1' * D * BL1 + G11 * eye(2), BL1' * D * BL2 + G12 * eye(2), BL1' * D * BL3 + G13 * eye(2), BL1' * D * BL4 + G14 * eye(2);
             BL2' * D * BL1 + G21 * eye(2), BL2' * D * BL2 + G22 * eye(2), BL2' * D * BL3 + G23 * eye(2), BL2' * D * BL4 + G24 * eye(2);
             BL3' * D * BL1 + G31 * eye(2), BL3' * D * BL2 + G32 * eye(2), BL3' * D * BL3 + G33 * eye(2), BL3' * D * BL4 + G34 * eye(2);
             BL4' * D * BL1 + G41 * eye(2), BL4' * D * BL2 + G42 * eye(2), BL4' * D * BL3 + G43 * eye(2), BL4' * D * BL4 + G44 * eye(2)];
    
    re_qp = [BL1' * Se;
             BL2' * Se;
             BL3' * Se;
             BL4' * Se];
         
    kte = kte + kte_qp * w(iqp) * det(JofXiEta);
    re = re + re_qp * w(iqp) * det(JofXiEta);
    fe = fe + (N' * BodyForce(XY)) * w(iqp) * det(JofXiEta);
end