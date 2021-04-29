% LoadStep()
% At each loadstep increment, this function calculates the change in 
% displacement and determines whether convergence has been reached
% NOTE: FUNCTION CALLS DO NOT INCLUDE INPUTS/OUTPUTS IN FINAL FORM
function [] = LoadStep()

% Load initial displacement vector
vk = v_nought;

% Loop through each loadstep
L = linspace(0,maxLoad,numLoadSteps + 1); 
for k = 1:1:size(L,1)
    % Initialize parameters for each loadstep
    loadstep = L(k); vi = vk; convcon = 0; iterations = 0; 
    
    % Find G, KT at first iteration using GlobalSystem.m
    [G,KT] = GlobalSystem(vi,loadstep);
    
    % Determine if norm(G) at first iteration is zero
    normGi = norm(G);
    if normGi == 0
        % Iterate loadstep
        convcon = 1;
    else
        iterations = iterations + 1;
    end
    
    % Until the convergence condition is satisfied, iterate the
    % displacement vector
    while (convcon == 0) && (iterations < 15)
        % Find iterated displacement v_i+1 using Soln.m
        deltav = Soln(G,KT); 
        vipo = vi + deltav;     
        
        % Find G, KT at iteration i+1 using GlobalSystem.m
        [G,KT] = GlobalSystem(vipo,loadstep);
        
        % Test for convergence: If norm(G_i+1) <= (10^-4)*norm(G_i), set 
        % v_k+1 = v_i+1 and move on to next loadstep (k+1)
        normGipo = norm(G);  
        if normGipo <= (normGi * (10^-4))
            vk = vipo; convcon = 1; 
        end

        % Update values for next iteration
        vi = vipo; normGi = normGipo;
        iterations = iterations + 1; 
    end
    
    % End of current loadstep, package iterated results
end