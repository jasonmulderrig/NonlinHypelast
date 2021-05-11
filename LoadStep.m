function [solverStruct,globalSystem] = LoadStep(meshStruct,boundStruct,solverStruct,globalSystem)
% [solverStruct,globalSystem] = LoadStep(meshStruct,boundStruct,solverStruct,globalSystem);
% Use a Newton-Raphson iterative solution scheme to solve for the global 
% displacement vector. At each increment, this function calculates the 
% change in displacement and determines if convergence has been reached
% NOTE: FUNCTION CALLS DO NOT INCLUDE INPUTS/OUTPUTS IN FINAL FORM

% last update: 30 Apr 2021 C. Bonneville; J. Mulderrig; S. Srivatsa 

% Initialization before looping through each load increment
numIncrements=solverStruct.numIncrements;
maxIterations=solverStruct.maxIterations;

iterations2conv = NaN(numIncrements,1);
convergenceBehavior = NaN(numIncrements,maxIterations+1);

d_i = globalSystem.d; % initial global displacement vector

tol = solverStruct.tol; % absolute convergence tolerance value

tic
% Loop through each load increment
for k = 1:numIncrements
    toc
    disp(sprintf('\b   (Nonlinear Solution of Load Step %d)',k))
    tic
    % Initialization for each load increment
    boundStruct.SurfEss = boundStruct.SurfEssIncrements{k};
    boundStruct.SurfNat = boundStruct.SurfNatIncrements{k};
    iterations = 0;
    
    % Apply the incremental essential boundary conditions at the initial
    % zeroth iteration. Note that the incremental essential boundary
    % conditions vary only if non-homogeneous essential boundary conditions
    % are applied.
    [d_i, boundStruct] = ApplyAllEssBCs(d_i,boundStruct);
    freeDOF = boundStruct.freeDOF;
    
    % Calculate G and K_T at the initial zeroth iteration
    [G_i,K_T_i] = GlobalSystemCalcn(d_i,meshStruct,boundStruct);
    
    % Determine if G is the zero vector at the initial zeroth iteration,
    % which will analytically occur at the first load increment (accounting
    % for zero applied force) with the initial zero global displacement
    % vector
    norm_G_i_F = norm(G_i(freeDOF));
    convergenceBehavior(k,iterations+1) = norm_G_i_F;
    if norm_G_i_F < tol
        iterations2conv(k) = iterations;
        continue;
    end
    
    % Until convergence is reached, iteratively solve for the converged
    % displacement vector
    while true
        iterations = iterations+1;
        % Find iterated displacement d_i_plus_1 using Soln.m
        delta_d_i_plus_1 = SolnNL(G_i,K_T_i,boundStruct);
        d_i_plus_1 = d_i + delta_d_i_plus_1;   
        
        % Find G, KT at iteration i+1 using GlobalSystem.m
        [G_i_plus_1,K_T_i_plus_1] = GlobalSystemCalcn(d_i_plus_1,meshStruct,boundStruct);
        
        norm_G_i_plus_1_F = norm(G_i_plus_1(freeDOF));
        
        % Test for convergence:
        if norm_G_i_plus_1_F <= tol % convergence is achieved
            d_i = d_i_plus_1; % update global displacement vector
            iterations2conv(k) = iterations;
            convergenceBehavior(k,iterations+1) = norm_G_i_plus_1_F;
            break;
        end
        % convergence is not achieved AND the maximum number of iterations
        % has not been reached
        if (( norm_G_i_plus_1_F > tol || isnan(norm_G_i_plus_1_F) ) && iterations < maxIterations)
            d_i = d_i_plus_1; % update global displacement vector
            G_i = G_i_plus_1; % update global residual vector
            K_T_i = K_T_i_plus_1; % update global tangent stiffness matrix
            convergenceBehavior(k,iterations+1) = norm_G_i_plus_1_F;
            continue;
        end
        % convergence is not achieved AND the maximum number of iterations
        % has been reached
        if (( norm_G_i_plus_1_F > tol || isnan(norm_G_i_plus_1_F) ) && iterations == maxIterations)
            d_i = d_i_plus_1; % update global displacement vector
            % indicate that convergence was not reached at the maximum
            % iteration number
            iterations2conv(k) = iterations+1;
            convergenceBehavior(k,iterations+1) = norm_G_i_plus_1_F;
            break;
        end
    end
end
globalSystem.d=d_i; % Store final converged global displacement vector
% Store vector of iterations needed to achieve convergence at each load
% increment
solverStruct.iterations2conv=iterations2conv;
% Store matrix containing the values of the L2 norm of the global residual
% vector for each Newton-Raphson iteration at each incremental load step
solverStruct.convergenceBehavior=convergenceBehavior;