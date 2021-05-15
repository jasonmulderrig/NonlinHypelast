function [d, boundStruct] = ApplyAllEssBCs(d,boundStruct)
% [d,boundStruct] = ApplyAllEssBCs(d,boundStruct)
% Apply all essential boundary conditions to the global displacement vector
% by caling the ApplyEssBC.m method. Distinguish and store which degrees of 
% freedom are and are not associated with the essential boundary conditions

% last update: 16 May 2021 C. Bonneville; J. Mulderrig; S. Srivatsa 

essDOF = []; % essDOF stores all the essential DOF
essVals = [];  % essVals stores all the corresponding applied essential values

for i = 1:size(boundStruct.SurfEss,1)   % Loop over all essential boundaries 

    [cessDOF, cessVals]=ApplyEssBC(i,boundStruct);% assign essential BCs
    essDOF =  [essDOF; cessDOF];
    essVals = [essVals; cessVals];

end

% Because the corner nodes are at two boundaries, we need to find the
% unique numbers of debc (you dont want to apply the same EBC twice!)
[essDOF,m,~] = unique(essDOF);
essVals = essVals(m);

% Store the DOF associated with the essential BCs
boundStruct.essDOF = essDOF;
d(essDOF) = essVals;

% Store the DOF that are not associated with the essential BCs
numEq=length(d);
freeDOF=setdiff(1:numEq,essDOF)';
boundStruct.freeDOF = freeDOF;


