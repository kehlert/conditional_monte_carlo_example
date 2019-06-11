function [X, freq, total_rvs_used] = cmc(T, x0, N, test_N)  
    % cmc Apply conditional Monte Carlo to the Lotka-Volterra model.
    %   example: [X, freq] = cmc(4, [200, 100], 10^3, 10)
    %
    %   T = float
    %   x0 = row vector
    %   N = integer
    %   test_N = integer
    %
    %   The optimal m and h are estimated based on test_N paths. Then
    %   conditional Monte Carlo is run with N independent 'branches'. In
    %   other words, we use the optimal m and h to obtains N*m observations
    %   of the state of the Lotka-Volterra model at time t.
    %
    %   Returns the states we observed in X, and a matrix of observed state
    %   frequencies in freq. Also returns the total number of random 
    %   variates used in the simulations.
    %
    %   NOTE: in order to run the simulations faster, we can use Matlab's
    %   coder functionality, which turns Matlab code into compiled C code.
    %   In order to enable this, run gen_c_code.m first.
   
    addpath('./permn')
    
    mex_detected = true;
    if ~exist('codegen', 'dir')
        fprintf(['WARNING: Did not detect the ./codegen directory. To ' ...
                 'make the simulations faster, run gen_c_code.m.\n']);
            
        mex_detected = false;
    end
    
    %sample intensities in order to estimate P(X_11 = X_12)
    fprintf('Sampling intensity paths...');

    lambdas_paths = cell(test_N,1);
    for n = 1:test_N
        [~, new_lambdas_path, new_lambdas_times] = sim_lotka_volterra(T, x0, true);
        lambdas_paths{n} = timeseries(new_lambdas_path, new_lambdas_times);
        lambdas_paths{n} = setinterpmethod(lambdas_paths{n}, 'linear');
        lambdas_paths{n}.Name = num2str(n);
    end

    %compute the average total intensity path
    new_times = linspace(0, T, 10001);
    number_of_reactions = size(lambdas_paths{1}.Data,2);
    resampled_lambda_paths = zeros(length(new_times), number_of_reactions, test_N);

    for n=1:test_N
        lambdas_paths{n} = resample(lambdas_paths{n}, new_times);
        resampled_lambda_paths(:,:,n) = lambdas_paths{n}.Data;
    end

    avg_total_intensity_path = timeseries(mean(sum(resampled_lambda_paths,2), 3), new_times);
    avg_total_intensity_path = setinterpmethod(avg_total_intensity_path, 'linear');
    avg_total_intensity_path.Name = 'average total intensity path';
    fprintf('done.\n');

    %find the optimal m and h based on avg_total_intensity_path
    first_term = @(m,h) cost_function(avg_total_intensity_path, T, m, h) ./ m;

    %Matlab doesn't support bound-constrained derivative-free optimization,
    %so I just came up with a quick hack to deal with that.
    %I essentially just set the value of the objective to 10^20 (which is very large
    %compared to typical values) if we are outside the bounds. It's not
    %pretty, but we are just trying to demonstrate the idea, not provide
    %the best possible implementation.
    obj_func = @(pair) (first_term(pair(1),pair(2)) * (1 + (pair(1)-1) * get_joint_prob(pair(2), lambdas_paths))) * (pair(1) >= 1) * (pair(2) >= 0) * (pair(2) <= T) + 10^20 * ((pair(1) < 1) + (pair(2) < 0) + (pair(2) > T));
    
    opt_m_h = fminsearch(obj_func, [100, T/10]);
    fprintf('Optimal (m,h): (%d, %f)\n', opt_m_h(1), opt_m_h(2));
    m = ceil(opt_m_h(1));
    h = opt_m_h(2);
    
    fprintf('Sampling X_{ij} for i=1:N and j=1:m...');
    total_rvs_used = 0;
    X = zeros(N * m, length(x0));

    for i = 1:N
        if mod(i, ceil(N/10)) == 0
            fprintf('\nsimulations run: %d', i);
        end
        
        if mex_detected
            [x, ~, ~, rvs_used] = sim_lotka_volterra_mex(T-h, x0, false);
        else
            [x, ~, ~, rvs_used] = sim_lotka_volterra(T-h, x0, false);
        end
        
        total_rvs_used = total_rvs_used + rvs_used;
        
        for j = 1:m
            if mex_detected
                [X((i-1) * m + j, :), ~, ~, rvs_used] = sim_lotka_volterra_mex(h, x, false);
            else
                [X((i-1) * m + j, :), ~, ~, rvs_used] = sim_lotka_volterra(h, x, false);
            end
            
            total_rvs_used = total_rvs_used + rvs_used;
        end
    end
    
    dist_ub = [800, 800];
    freq = zeros(dist_ub+1);
    for i=1:length(X)
        %add 1 in case X contains a zero
        if all(X(i,:) <= dist_ub)
            freq(X(i,1)+1,X(i,2)+1) = freq(X(i,1)+1,X(i,2)+1) + 1;
        end
    end
    
    freq = freq / (N * m);

    fprintf('\ndone.\n');
end
