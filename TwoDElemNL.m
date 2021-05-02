function [K_T_e, fe, re] = TwoDElemNL(meshStruct, elmID, qp, w, d)
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
ConstitutiveLaw=meshStruct.ConstitutiveLaw;

D = [lambda + 2 * mu, lambda, 0;
     lambda, lambda + 2 * mu, 0;
     0, 0, mu]; % St Venant formulation

% get nodal coordinates for this element
xvec = nCoords(elCon(elmID,:), 1); % these are column vectors of the
yvec = nCoords(elCon(elmID,:), 2); % physical coordinates of the
                                % nodes of this element

% nodal x-coordinate in current configuration
ux = xvec + d(elCon(elmID, :) * 2 - 1); 
% nodal y-coordinate in current configuration
uy = yvec + d(elCon(elmID, :) * 2);

K_T_e = zeros(numDOF * nnpe);  % initialize elemental stiffness matrix
G_e = zeros(numDOF * nnpe, 1); % initialize elemental residual vector 
% fe = zeros(numDOF * nnpe, 1); % initialize elemental force vector
% re = zeros(numDOF * nnpe, 1); % initialize elemental residual vector

for iqp = 1 : size(qp, 1) % loop over quadrature points
    % for each quadrature point ...
    NofXiEta = Nmatrix(qp(iqp, :), nnpe); % row vector of shape functions
                                      % at this QP in parent domain
    dNdXiEta = dNmatrix(qp(iqp, :),nnpe);% 2xnnpe array of shape func derivs
                                      % at this QP in parent domain
    JofXiEta = dNdXiEta * [xvec, yvec];    % jacobian at this QP (2x2 matrix)
    dNdXY = inv(JofXiEta) * dNdXiEta;     % 2xnnpe array of dNdX at ths QP
    
    XY = NofXiEta * [xvec, yvec]; % physical coordinate of this QP (1x2)
    
    % now create the N matrix
    N = [];
    for np=1:nnpe
        N = [N,[NofXiEta(np), 0; 0, NofXiEta(np)]];
    end
    
    % -------- Deformation gradient F_e (Slide 13/18 - Eq 24) ----------- %
    
    F_e = zeros(2);
    for np=1:nnpe
        F_e_np = [ux(np) * dNdXY(1, np), ux(np) * dNdXY(2, np);
                  uy(np) * dNdXY(1, np), uy(np) * dNdXY(2, np)];
        F_e = F_e + F_e_np; 
    end
    
    % ------------- Right Cauchy-Green strain tensor C_e ---------------- %
    
    C_e = Fe' * Fe;
    
    % --------- Euler Lagrange Tensor Ee (Slide 13/18 - Eq 25) ---------- %
    
    E_e = 0.5 * (C_e - eye(2));
    
    % ------------------- Material stiffness matrix D ------------------- %
    
    switch ConstitutiveLaw
        case 'compressibleNeoHookean'
            % solve for D here for the element and qp, which is dependent 
            % on C_e and det(F_e)
        case 'StVenant'
            D=meshStruct.Material.D;
    end
    
    % ----------------- B_LI Matrix (Slide 14/18 - Eq 26) --------------- %
    
    B_LI_e = cell(nnpe,1);
    for np=1:nnpe
        B_LI_e(np) = [F_e(1, 1) * dNdXY(1, np), F_e(2, 1) * dNdXY(1, np);
                      F_e(1, 2) * dNdXY(2, np), F_e(2, 2) * dNdXY(2, np);
                      F_e(1, 1) * dNdXY(2, np) + F_e(1, 2) * dNdXY(1, np),  Fe(2, 1) * dNdXY(2, np) + Fe(2, 2) * dNdXY(1, np)];
    end
       
    % --------------- Second PK Stress (Slide 12/18 - Eq 23) ------------ %
    
    S_e = D * [E_e(1, 1); E_e(2, 2); 2 * E_e(1, 2)]; % Voigt notation
    
    % ------------------ G Matrix (Slide 16/18 - Eq 31) ---------------- %
    
    S_e_matrix = [S_e(1), S_e(3); S_e(3), S_e(2)]; % Matrix notation
    
    G_e = zeros(nnpe);
    for I=1:nnpe
        for K=1:nnpe
            G_e(I,K) = dNdXY(:, I)' * S_e_matrix * dNdXY(:, K);
        end
    end
    
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental body force vector
    
    kte_qp = [BL1' * D * BL1 + G11 * eye(2), BL1' * D * BL2 + G12 * eye(2), BL1' * D * BL3 + G13 * eye(2), BL1' * D * BL4 + G14 * eye(2);
             BL2' * D * BL1 + G21 * eye(2), BL2' * D * BL2 + G22 * eye(2), BL2' * D * BL3 + G23 * eye(2), BL2' * D * BL4 + G24 * eye(2);
             BL3' * D * BL1 + G31 * eye(2), BL3' * D * BL2 + G32 * eye(2), BL3' * D * BL3 + G33 * eye(2), BL3' * D * BL4 + G34 * eye(2);
             BL4' * D * BL1 + G41 * eye(2), BL4' * D * BL2 + G42 * eye(2), BL4' * D * BL3 + G43 * eye(2), BL4' * D * BL4 + G44 * eye(2)];
    
    re_qp = [BL1' * S_e;
             BL2' * S_e;
             BL3' * S_e;
             BL4' * S_e];
         
    K_T_e = K_T_e + kte_qp * w(iqp) * det(JofXiEta);
    re = re + re_qp * w(iqp) * det(JofXiEta);
    fe = fe + (N' * BodyForce(XY)) * w(iqp) * det(JofXiEta);
end