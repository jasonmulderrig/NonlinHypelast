% NonlinHypelast
%
% This code solves 2D nonlinear hyperelasticity problems using 3-node
% or 4-node elements. It is heavily based on N. Zabaras's code, and
% modified by Christophe Bonneville, Jason Mulderrig, and Srikar Srivatsa
% in the spring of 2021.
%
% The list of functions which must be edited for each problem is:
% MeshGeometry      : define the domain, element type, set plotting 
%                     options, plot the mesh, etc.
% InputData         : define natural and essential BC values, define
%                     material stiffness matrix D.
% BodyForce         : define x and y components of body force.
% 
% All of the structs used, including their fields and the
% functions in which they are first defined, are listed here:
% 
% meshStruct : defined in MeshGeometry
%     meshStruct.elCon     : connectivity array
%     meshStruct.gatherMat : gather/scatter matrix
%     meshStruct.nCoords   : nodal coordinates
%     meshStruct.nnpe      : number of nodes per element
%     meshStruct.nsd       : number of spatial dimensions
%     meshStruct.numDOF    : number of degrees of freedom per node
%     meshStruct.numEls    : number of elements
%     meshStruct.numEq     : number of equations in global system
%     meshStruct.numNodes  : number of global nodes
%     meshStruct.Material  : struct with material properties defined in
%                            InputData
%         meshStruct.Material.E : Young's modulus  
%         meshStruct.Material.nu: Poisson's ratio  
%         meshStruct.Material.D : constitutive relation matrix 
%                                 such that stress=D*strain  
%
% boundStruct : defined in BoxGrid_2D (or loadFromGridFile) and InputData
%     boundStruct.elements
%         boundStruct.elements.Elems            : element numbers
%                                                    on each boundary
%         boundStruct.elements.SurfaceIndicator : flag for which element
%                                                    edge is on boundary
%     boundStruct.nodes
%         boundStruct.nodes.Nodes               : global node number on
%                                                    each boundary
%     boundStruct.SurfNat                       : prescribed natural BC
%     boundStruct.SurfEssV                      : prescribed essential BC
% 
% PlotInstructions : defined in InputData
%     PlotInstructions.plot_mesh     : flag for plotting mesh
%     PlotInstructions.plot_node     : flag for plotting node numbers
%     PlotInstructions.plot_boundary : flag for plotting boundary numbers
%     PlotInstructions.plot_contour  : flag for plotting nodal solution
%     PlotInstructions.plot_deformed : flag for plotting deformed mesh
%     PlotInstructions.plot_fringes  : flag for plotting fringe pattern

% last update: 16 May 2021 C. Bonneville; J. Mulderrig; S. Srivatsa 

% preliminary necessities
home; clear; close all; 
tic % start a timer, just for fun

% set up the mesh and define the boundaries
[meshStruct,boundStruct,PlotInstructions]=MeshGeometry;

% choose plotting options and define BC type and material properties
[meshStruct,boundStruct,solverStruct,globalSystem]=InputData(meshStruct,boundStruct);
toc
disp(sprintf('\b   (Preprocessing)')) % output the time for meshing
tic
% Solve the global displacements using a Newton-Raphson iterative solver
globalSystem = LoadStep(meshStruct,boundStruct,solverStruct,globalSystem);
toc
disp(sprintf('\b   (Nonlinear Solution)')) % output the time for the solution
tic

% Post-process for plots, flux, etc.
PostProcessor(PlotInstructions,meshStruct,globalSystem);
toc
disp(sprintf('\b   (Postprocessing)')) % output time for post-processing
