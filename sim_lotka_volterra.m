function [x, lambdas_path, lambdas_times, rvs_used] = sim_lotka_volterra(T, x0, get_lambdas_path)
    % sim_lotka_volterra Simulate a single path of the Lotka-Volterra model.
    %   example: [x, lambdas_path, lambdas_times] = sim_lotka_volterra(4, 200, 100], false)
    %
    %   T = float
    %   x0 = row vector
    %   get_lambdas_path = true/false
    %
    %   Simulates the Lotka-Volterra model and returns the state x at time
    %   T. If get_lambdas_path is true, then the function also returns the
    %   intensities of the reactions after every event in lambdas_path, and the times of
    %   the events in lambdas_times.

	theta = [2, 0.01, 2];

    x = x0;
	t = 0;

	lambdas = [theta(1)*x(1), theta(2)*x(1)*x(2), theta(3)*x(2)];

    coder.varsize('lambdas_times', [10^8, 1], [1 0]);
	lambdas_times = zeros(1);
    
    coder.varsize('lambdas_path', [10^8, 3], [1 0]);
    lambdas_path = zeros(1, length(lambdas));
	lambdas_path(1,:) = lambdas;
    
    rvs_used = 0;

    while true
	    lambda0 = sum(lambdas);
	    dt = exprnd(1/lambda0);
	    
	    if t+dt > T
            break;
	    end
	    
	    t = t + dt;
	    
	    u = rand;
	    if u < lambdas(1) / lambda0
            x(1) = x(1) + 1;
	    elseif u < (lambdas(1)+lambdas(2)) / lambda0
            x(1) = x(1) - 1;
            x(2) = x(2) + 1;
        else
            x(2) = x(2) - 1;
	    end
	    
	    lambdas = [theta(1)*x(1), theta(2)*x(1)*x(2), theta(3)*x(2)];
	    
        if get_lambdas_path
            lambdas_path(size(lambdas_path,1) + 1,:) = lambdas;
            lambdas_times(length(lambdas_times) + 1) = t;
        end
        
        rvs_used = rvs_used + 2;
    end
    
    if get_lambdas_path
        lambdas_path(size(lambdas_path,1) + 1, :) = lambdas;
        lambdas_times(length(lambdas_times) + 1) = T;
    end
end
