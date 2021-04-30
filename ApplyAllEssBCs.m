% ApplyAllEssBCs
% This function applies essential BCs (calling subfunction ApplyEssBC.m). 

% last update: 29 Apr 2021 J. Mulderrig  

function boundStruct = ApplyAllEssBCs(boundStruct)

% Apply essential BC(s)

essDOF = []; % essDOF stores all the essential DOF
essVals = [];  % essVals stores all the corresponding applied essential values

for i = 1:size(boundStruct.SurfEssV,1)   % Loop over all essential boundaries 

    [cessDOF, cessVals]=ApplyEssBC(i,boundStruct);% assign essential BCs
    essDOF =  [essDOF; cessDOF];
    essVals = [essVals; cessVals];

end

% Because the corner nodes are at two boundaries, we need to find the
% unique numbers of debc (you dont want to apply the same EBC twice!)
[essDOF,m,~] = unique(essDOF);
essVals = essVals(m);

boundStruct.essDOF = essDOF;
boundStruct.ebcVals = essVals;
