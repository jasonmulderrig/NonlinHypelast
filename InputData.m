function [meshStruct,boundStruct,solverStruct,globalSystem]=InputData(meshStruct,boundStruct)
% [meshStruct,boundStruct,solverStruct,globalSystem]=InputData(meshStruct,boundStruct);
% Define the essential and natural BCs, and define the material stiffness
% matrix D

% last update: 30 Apr 2021 C. Bonneville; J. Mulderrig; S. Srivatsa 

% For complex geometries, you need to firstly plot the mesh (set both 
% PlotInstructions.plot_mesh and PlotInstructions.plot_boundary to be
% 'yes') to see how the number of boudaries are defined.

% Define the essential BCs
% boundStruct.SurfEss = [];
boundStruct.SurfEss = [2 1 0
                        2 2 0
                        4 1 0.5]; % e.g. [4 2 20] means all nodes on surface # 4,
%                                % degree of freedom #2 (y direction), has a value of 20.

% Define the natural BCs
% The natural boundary condition is defined in tangential and normal
% direction (rather than global x and y direction),outer normal is
% positive.
boundStruct.SurfNat = [3 0 1e11]; % e.g. [3 10 -10] means surface # 3 has 
                                 % a constantly distributed tangential traction
                                 % 10 and normal traction (pointing in) 10.

% Define material properties
Lambda = 1; % first normalized Lame constant
mu     = 10; % second normalized Lame constant
DeformationState = 'PlaneStrain'; % only in plane strain deformation
ConstitutiveLaw = 'StVenant';
switch DeformationState
    case 'PlaneStrain' % do nothing; continue on
    otherwise
        error('Is this two-dimensional plane strain deformation?');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ConstitutiveLaw
    case 'StVenant'
        D=[Lambda+2*mu Lambda 0 ; Lambda Lambda+2*mu 0 ; 0 0 mu];
        meshStruct.Material.D=D;
    case 'compressibleNeoHookean' % do nothing; continue on
    otherwise
        error('Is this the St Venant or the compressible Neo-Hookean constitutive law?');
end

meshStruct.Material.Lambda=Lambda;
meshStruct.Material.mu=mu;
meshStruct.DeformationState=DeformationState;
meshStruct.ConstitutiveLaw=ConstitutiveLaw;

% initialize the global displacement vector to the zero vector
numEq=MeshStruct.numEq;
d=zeros(numEq,1);
globalSystem.d = d;

numIncrements = 41; ement

% Define the displacement increments for all the essemtial boundary conditions
SurfEssIncrements = cell(numIncrements,1);

for k = 1:numIncrements % Loop over all displacement increments
    SurfEssIncrement = boundStruct.SurfEss;
    for i = 1:size(SurfEssIncrement,1) % Loop over all essential boundaries
        fixedDisplacementVal = SurfEssIncrement(i,3);
        fixedDisplacementIncr = (k-1)/(numIncrements-1)*fixedDisplacementVal;
        SurfEssIncrement(i,3) = fixedDisplacementIncr;
    end
    SurfEssIncrements(k) = SurfEssIncrement;
end
boundStruct.SurfEssIncrements = SurfEssIncrements;

% Define the load increments for all the natural boundary conditions
SurfNatIncrements = cell(numIncrements,1);

for k = 1:numIncrements % Loop over all load increments
    SurfNatIncrement = boundStruct.SurfNat;
    for i = 1:size(SurfNatIncrement,1) % Loop over all natural boundaries
        tangentialTractionVal = SurfNatIncrement(i,2);
        normalTractionVal = SurfNatIncrement(i,3);
        tangentialTractionIncr = (k-1)/(numIncrements-1)*tangentialTractionVal;
        normalTractionIncr = (k-1)/(numIncrements-1)*normalTractionVal;
        SurfNatIncrement(i,2) = tangentialTractionIncr;
        SurfNatIncrement(i,3) = normalTractionIncr;
    end
    SurfNatIncrements(k) = SurfNatIncrement;
end
boundStruct.SurfNatIncrements = SurfNatIncrements;

% Define constants used for the nonlinear Newton-Raphson solution scheme
solverStruct.numIncrements=numIncrements; % number of 
% maximum number of Newton-Raphson iterations permitted
solverStruct.maxIterations = 15; 
% used to define the relative convergence tolerance
solverStruct.epsilon = 1e-4; 


