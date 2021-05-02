function globalSystem = LoadStep(meshStruct,boundStruct,solverStruct,globalSystem)
% globalSystem = LoadStep(meshStruct,boundStruct,solverStruct,globalSystem);
% Use a Newton-Raphson iterative solution scheme to solve for the global 
% displacement vector. At each increment, this function calculates the 
% change in displacement and determines if convergence has been reached
% NOTE: FUNCTION CALLS DO NOT INCLUDE INPUTS/OUTPUTS IN FINAL FORM

% last update: 30 Apr 2021 C. Bonneville; J. Mulderrig; S. Srivatsa 

% Initialization before looping through each load increment
numIncrements=solverStruct.numIncrements;
maxIterations=solverStruct.maxIterations;

iterations2conv = NaN(numIncrements,1);
d_i = globalSystem.d; % initial global displacement vector

% Loop through each load increment
for k = 1:numIncrements
    % Initialization for each load increment
    boundStruct.SurfEss = boundStruct.SurfEssIncrs{k};
    boundStruct.SurfNat = boundStruct.SurfNatIncrs{k};
    iterations = 0;
    
    % Apply the incremental essential boundary conditions at the initial
    % zeroth iteration. Note that the incremental essential boundary
    % conditions vary only if non-homogeneous essential boundary conditions
    % are applied.
    boundStruct = ApplyAllEssBCs(boundStruct);
    
    % Calculate G and K_T at the initial zeroth iteration
    [G_i,K_T_i] = GlobalSystemCalcn(d_i,meshStruct,boundStruct);
    
    % Determine if G is the zero vector at the initial zeroth iteration,
    % which will analytically occur at the first load increment (accounting
    % for zero applied force) with the initial zero global displacement
    % vector
    if norm(G_i) == 0
        iterations2conv(k) = iterations;
        continue;
    else
        epsilon = solverStruct.epsilon;
        tol = epsilon*norm(G_i);
    end
    
    % Until convergence is reached, iteratively solve for the converged
    % displacement vector
    while true
        iterations = iterations+1; 
        % Find iterated displacement d_i_plus_1 using Soln.m
        delta_d_i_plus_1 = Soln(G_i,K_T_i); 
        d_i_plus_1 = d_i + delta_d_i_plus_1;     
        
        % Find G, KT at iteration i+1 using GlobalSystem.m
        [G_i_plus_1,K_T_i_plus_1] = GlobalSystemCalcn(d_i_plus_1,meshStruct,boundStruct);
        
        % Test for convergence:
        if norm(G_i_plus_1) <= tol % convergence is achieved
            d_i = d_i_plus_1; % update global displacement vector
            iterations2conv(k) = iterations;
            break;
        end
        % convergence is not achieved AND the maximum number of iterations
        % has not been reached
        if ( norm(G_i_plus_1) <= tol || isnan(norm(G_i_plus_1)) ) && iterations < maxIterations
            d_i = d_i_plus_1; % update global displacement vector
            G_i = G_i_plus_1; % update global residual vector
            K_T_i = K_T_i_plus_1; % update global tangent stiffness matrix
            continue;
        end
        % convergence is not achieved AND the maximum number of iterations
        % has been reached
        if ( norm(G_i_plus_1) <= tol || isnan(norm(G_i_plus_1)) ) && iterations == maxIterations
            d_i = d_i_plus_1; % update global displacement vector
            % indicate that convergence was not reached at the maximum
            % iteration number
            iterations2conv(k) = iterations+1;
            break;
        end
    end
end
globalSystem.d=d_i; % Store final converged global displacement vector
% Store vector of iterations needed to achieve convergence at each load
% increment
solverStruct.iterations2conv=iterations2conv;