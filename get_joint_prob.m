function [joint_prob] = get_joint_prob(h, lambdas_paths)
    joint_prob = 0;
    nullspace = [1; 1; 1];
    dim_nullspace = 1;
    
    T = lambdas_paths{1}.Time(end);
        
    %compute the integral of average lambdas
    lambda_integrals = zeros(length(lambdas_paths), size(lambdas_paths{1}.Data, 2));

    for n = 1:length(lambdas_paths)
        path_subset = getsampleusingtime(lambdas_paths{n}, T-h, T);

        if length(path_subset.Time) == 1
            lambda_integrals(n,:) = h * path_subset.Data;
        else
            lambda_integrals(n,:) = trapz(path_subset.Time, path_subset.Data);
        end
    end

    max_k = 4;

    if isempty(nullspace) || max_k == 0
        joint_prob = mean(prod(max(besseli(0, 2 * lambda_integrals, 1), 0), 2));
    else
        coeffs = permn(-max_k:max_k, dim_nullspace);

        for i=1:length(coeffs)
            k = nullspace * coeffs(i,:)';
            rep_k = repmat(k', length(lambdas_paths), 1);
            skellam_prob_vec = prod(max(besseli(rep_k, 2 * lambda_integrals, 1), 0), 2);

            %besseli() gives extremely wrong results for large
            %parameter values, so we just set them to zero, since they
            %basically should be essentially zero anyway
            zero_indices = sum(abs(rep_k) > 5 * lambda_integrals, 2) > 0;
            skellam_prob_vec(zero_indices) = 0;

            skellam_prob = mean(skellam_prob_vec);
            joint_prob = joint_prob + skellam_prob;
        end
    end
end

