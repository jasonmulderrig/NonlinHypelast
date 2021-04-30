function [meshStruct,boundStruct,solverStruct,globalSystem]=InputData(meshStruct,boundStruct)
% [meshStruct,boundStruct,solverStruct,globalSystem]=InputData(meshStruct,boundStruct);
% Define the essential and natural BCs, and define the material stiffness
% matrix D

% last update: 15 Nov 2015 H. Ritz; Y. Xu

% For complex geometries, you need to firstly plot the mesh (set both 
% PlotInstructions.plot_mesh and PlotInstructions.plot_boundary to be
% 'yes') to see how the number of boudaries are defined.

% Define the essential BCs
% boundStruct.SurfEssV = [];
boundStruct.SurfEssV = [2 1 0
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

numIncrs = 41;

% Define the load increments for all the natural boundary conditions
numSurfNatIncrs = numIncrs;
SurfNatIncrs = cell(numSurfNatIncrs,1);

for h = 1:numSurfNatIncrs % Loop over all load increments
    SurfNatIncr = boundStruct.SurfNat;
    for i = 1:size(SurfNatIncr,1) % Loop over all natural boundaries
        tangentialTractionVal = SurfNatIncr(i,2);
        normalTractionVal = SurfNatIncr(i,3);
        tangentialTractionIncr = (h-1)/(numSurfNatIncrs-1)*tangentialTractionVal;
        normalTractionIncr = (h-1)/(numSurfNatIncrs-1)*normalTractionVal;
        SurfNatIncr(i,2) = tangentialTractionIncr;
        SurfNatIncr(i,3) = normalTractionIncr;
    end
    SurfNatIncrs(h) = SurfNatIncr;
end
boundStruct.SurfNatIncrs = SurfNatIncrs;

% Define the displacement increments for all the essemtial boundary conditions
numSurfEssIncrs = numIncrs;
SurfEssIncrs = cell(numSurfEssIncrs,1);

for h = 1:numSurfEssIncrs % Loop over all displacement increments
    SurfEssIncr = boundStruct.SurfEssV;
    for i = 1:size(SurfEssIncr,1) % Loop over all essential boundaries
        fixedDisplacementVal = SurfEssIncr(i,3);
        fixedDisplacementIncr = (h-1)/(numSurfEssIncrs-1)*fixedDisplacementVal;
        SurfEssIncr(i,3) = fixedDisplacementIncr;
    end
    SurfEssIncrs(h) = SurfEssIncr;
end
boundStruct.SurfEssIncrs = SurfEssIncrs;

% Define constants used for the nonlinear Newton-Raphson solution scheme
% maximum number of Newton-Raphson iterations permitted
solverStruct.maxIter = 15; 
% used to define the relative convergence tolerance
solverStruct.epsilon = 1e-4; 


