function [X, freq, total_rvs_used] = mc(T, x0, rvs_budget)  
    % mc Apply Monte Carlo to the Lotka-Volterra model.
    %   example: [X, freq] = mc(4, [200, 100], 10^6)
    %
    %   T = float
    %   x0 = row vector
    %   rvs_budget = integer
    %
    %   Runs Monte Carlo simulations until the number of random variates
    %   used by the simulations exceeds rvs_budget. The states observed are
    %   returned in X and the frequencies of the states are returned in
    %   freq.
    %
    %   NOTE: in order to run the simulations faster, we can use Matlab's
    %   coder functionality, which turns Matlab code into compiled C code.
    %   In order to enable this, run gen_c_code.m first.

    mex_detected = true;
    if ~exist('codegen', 'dir')
        fprintf(['WARNING: Did not detect the ./codegen directory. To ' ...
                 'make the simulations faster, run gen_c_code.m.\n']);
            
        mex_detected = false;
    end
    
    total_rvs_used = 0;
    [~, ~, ~, rvs_used] = sim_lotka_volterra_mex(T, x0, false);
    N_guess = ceil(rvs_budget / rvs_used);
    fprintf('Estimated number of simulations that need to run: %d', N_guess);
    X = zeros(N_guess, length(x0));
    X_index = 1;
    
    while total_rvs_used < rvs_budget
        if mod(X_index, ceil(N_guess/10)) == 0
            fprintf('\nsimulations run: %d', X_index);
        end
        
        if mex_detected
            [X(X_index, :), ~, ~, rvs_used] = sim_lotka_volterra_mex(T, x0, false);
        else
            [X(X_index, :), ~, ~, rvs_used] = sim_lotka_volterra(T, x0, false);
        end
            
        total_rvs_used = total_rvs_used + rvs_used;
        X_index = X_index + 1;
    end
    
    sims_run = X_index - 1;
    X = X(1:sims_run,:);
    
    dist_ub = [800, 800];
    freq = zeros(dist_ub+1);
    for i=1:length(X)
        %add 1 in case X contains a zero
        if all(X(i,:) <= dist_ub)
            freq(X(i,1)+1,X(i,2)+1) = freq(X(i,1)+1,X(i,2)+1) + 1;
        end
    end
    
    freq = freq / sims_run;

    fprintf('\ndone.\n');
end
