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
boundStruct.SurfEss = [4 1 0
                       4 2 0]; % e.g. [4 2 20] means all nodes on surface # 4,
%                                % degree of freedom #2 (y direction), has a value of 20.

% Define the natural BCs
% The natural boundary condition is defined in tangential and normal
% direction (rather than global x and y direction),outer normal is
% positive.
boundStruct.SurfNat = [2 0 1e-7]; % e.g. [3 10 -10] means surface # 3 has 
                                 % a constantly distributed tangential traction
                                 % 10 and normal traction (pointing in) 10.

% Define material properties
lambda = 1; % first normalized Lame constant
mu     = 10; % second normalized Lame constant
DeformationState = 'PlaneStrain'; % only in plane strain deformation
ConstitutiveLaw = 'compressibleNeoHookean';
switch DeformationState
    case 'PlaneStrain' % do nothing; continue on
    otherwise
        error('Is this two-dimensional plane strain deformation?');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ConstitutiveLaw
    case 'StVenant'
        D=[lambda+2*mu lambda 0 ; lambda lambda+2*mu 0 ; 0 0 mu];
        meshStruct.Material.D=D;
    case 'compressibleNeoHookean' % do nothing; continue on
    otherwise
        error('Is this the St Venant or the compressible Neo-Hookean constitutive law?');
end

meshStruct.Material.lambda=lambda;
meshStruct.Material.mu=mu;
meshStruct.DeformationState=DeformationState;
meshStruct.ConstitutiveLaw=ConstitutiveLaw;

% initialize the global displacement vector to the zero vector
numEq=meshStruct.numEq;
d=zeros(numEq,1);
globalSystem.d = d;

% number of increments for load and displacement permitted
numIncrements = 3; % 41; 

% Define the displacement increments for all the essemtial boundary conditions
SurfEssIncrements = cell(numIncrements,1);

for k = 1:numIncrements % Loop over all displacement increments
    SurfEssIncrement = boundStruct.SurfEss;
    for i = 1:size(SurfEssIncrement,1) % Loop over all essential boundaries
        fixedDisplacementVal = SurfEssIncrement(i,3);
        fixedDisplacementIncr = (k-1)/(numIncrements-1)*fixedDisplacementVal;
        SurfEssIncrement(i,3) = fixedDisplacementIncr;
    end
    SurfEssIncrements{k} = SurfEssIncrement;
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
    SurfNatIncrements{k} = SurfNatIncrement;
end
boundStruct.SurfNatIncrements = SurfNatIncrements;

% Define constants used for the nonlinear Newton-Raphson solution scheme
% number of increments for load and displacement permitted
solverStruct.numIncrements=numIncrements;
% maximum number of Newton-Raphson iterations permitted
solverStruct.maxIterations = 15; 
% used to define the constant convergence tolerance
solverStruct.tol = 1e-8; 


