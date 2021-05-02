function [P_e, R_e, K_T_e] = TwoDNLElem(d, meshStruct, elmID, qp, w)
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

% get nodal coordinates for this element
xvec = nCoords(elCon(elmID,:), 1); % these are column vectors of the
yvec = nCoords(elCon(elmID,:), 2); % physical coordinates of the
                                % nodes of this element

% nodal x-coordinate in current configuration
ux = xvec + d(elCon(elmID, :) * 2 - 1); 
% nodal y-coordinate in current configuration
uy = yvec + d(elCon(elmID, :) * 2);

P_e = zeros(numDOF * nnpe, 1); % initialize elemental external work vector
R_e = zeros(numDOF * nnpe, 1); % initialize elemental internal work vector
K_T_e = zeros(numDOF * nnpe);  % initialize elemental stiffness matrix

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
    
    C_e = F_e' * F_e;
    
    % --------- Euler Lagrange Tensor Ee (Slide 13/18 - Eq 25) ---------- %
    
    E_e = 0.5 * (C_e - eye(2));
    
    % ------------------- Material stiffness matrix D ------------------- %
    
    switch ConstitutiveLaw
        case 'compressibleNeoHookean'
            % solve for D here for the element and qp, which is dependent 
            % on C_e, det(F_e), and the Lame constants
            % 
            %
            %
        case 'StVenant'
            D=meshStruct.Material.D;
    end
    
    % ----------------- B_L Matrix (Slide 14/18 - Eq 26) ---------------- %
    
    B_L_e = cell(nnpe,1);
    for np=1:nnpe
        B_L_e(np) = [F_e(1, 1) * dNdXY(1, np), F_e(2, 1) * dNdXY(1, np);
                      F_e(1, 2) * dNdXY(2, np), F_e(2, 2) * dNdXY(2, np);
                      F_e(1, 1) * dNdXY(2, np) + F_e(1, 2) * dNdXY(1, np),  Fe(2, 1) * dNdXY(2, np) + Fe(2, 2) * dNdXY(1, np)];
    end
       
    % --------------- Second PK Stress (Slide 12/18 - Eq 23) ------------ %
    
    S_e = D * [E_e(1, 1); E_e(2, 2); 2 * E_e(1, 2)]; % Voigt notation
    
    % ------------------ G Matrix (Slide 16/18 - Eq 31) ---------------- %
    
    S_e_matrix = [S_e(1), S_e(3); S_e(3), S_e(2)]; % Matrix notation
    
    G_e_matrix = zeros(nnpe);
    for I=1:nnpe
        for K=1:nnpe
            G_e_matrix(I,K) = dNdXY(:, I)' * S_e_matrix * dNdXY(:, K);
        end
    end
    
    % add the weighted contributions of this QP to the elemental stiffness
    % matrix and elemental internal work vector
    
    % calculate elemental external work vector contribution of this QP from
    % the external applied body forces
    P_e_qp = N' * BodyForce(XY);
    
    % initialize elemental internal work vector contribution of this QP
    R_e_qp = zeros(numDOF * nnpe, 1);
    % initialize elemental stiffness matrix contribution of this QP
    K_T_e_qp = zeros(numDOF * nnpe);
    
    for I=1:nnpe
        i = 2 * I - 1; ii = 2 * I;
        B_LI_e = B_L_e{I};
        R_e_qp(i:ii) = B_LI_e' * S_e;
        
        for K=1:nnpe
            k = 2 * K - 1; kk = 2 * K; 
            B_LK_e = B_L_e{K};
            K_T_e_qp(i:ii,k:kk) = B_LI_e' * D * B_LK_e + G_e_matrix(I,K) * eye(2);
        end
    end
    
    
    P_e = P_e + P_e_qp * w(iqp) * det(JofXiEta);
    R_e = R_e + R_e_qp * w(iqp) * det(JofXiEta);
    K_T_e = K_T_e + K_T_e_qp * w(iqp) * det(JofXiEta);
end