% PostProcessor(PlotInstructions,meshStruct,globalSystem)
% Calculate stress, make contour plots, vector field plots, deformed mesh, etc.
%
% last update: 22 November 2017 H. Ritz
function PostProcessor(PlotInstructions,meshStruct,globalSystem,solverStruct)
d = globalSystem.d;
neq=meshStruct.numEq;
u = d(1:2:neq-1);
v = d(2:2:neq);

set(0,'defaultLineLineWidth',1.5)
set(0,'defaultTextFontSize',12)
set(0,'defaultAxesFontSize',12)
% set(0,'defaultAxesFontWeight','bold')


if strcmp(PlotInstructions.plot_deformed,'yes')
    PlotDeformedMesh('Deformed Configuration',meshStruct,d);
end
% Plot the contour distribution
if strcmp(PlotInstructions.plot_contour,'yes')
    patchPlot(u,meshStruct, 'FE solution: u displacement ');
    patchPlot(v,meshStruct, 'FE solution: v displacement ');
end

% Calculate the strain and stress
[~,sig] = calNodalStrainStress(d,meshStruct,globalSystem);

sigma_xx = sig(:,1);
sigma_yy = sig(:,2); 
sigma_xy = sig(:,3);

% Calculate the principal stress
sigma_1 = (sigma_xx + sigma_yy)/2 + sqrt(((sigma_xx + sigma_yy)/2).^2 + sigma_xy.^2);
sigma_2 = (sigma_xx + sigma_yy)/2 - sqrt(((sigma_xx - sigma_yy)/2).^2 + sigma_xy.^2);
tau_max = (sigma_1-sigma_2)/2;
sigma_vm=sqrt(sigma_xx.^2+sigma_yy.^2-sigma_xx.*sigma_yy+3*sigma_xy.^2);

% Plot stresses
if strcmp(PlotInstructions.plot_contour,'yes')

    patchPlot(sigma_xx,meshStruct, 'xx stress ');
    patchPlot(sigma_yy,meshStruct, 'yy stress ');
    patchPlot(sigma_xy,meshStruct, 'xy stress ');
    patchPlot(sigma_vm,meshStruct, 'von Mises stress ');
end

if strcmp(PlotInstructions.plot_fringes,'yes')
    PlotFringes(2*tau_max, 10,meshStruct);
end

% Plot Convergence Behavior
convBeh = solverStruct.convergenceBehavior;
iter = solverStruct.iterations2conv; 
numInc = solverStruct.numIncrements;
iter2inc = 0:1:iter(2); iterlastinc = 0:1:iter(numInc);
figure(9)
subplot(2,1,1)
plot(iter2inc,convBeh(2,1:(iter(2)+1)),'b-');
title('Convergence Behavior for Second Loadstep')
xlabel('Iteration Number')
ylabel('Norm of Free DOFs of Global Residual Vector')
set(gca,'YScale','log')
subplot(2,1,2)
plot(iterlastinc,convBeh(numInc,1:(iter(numInc)+1)),'b-');
title('Convergence Behavior for Last Loadstep')
xlabel('Iteration Number')
ylabel('Norm of Free DOFs of Global Residual Vector')
set(gca,'YScale','log')
