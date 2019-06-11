function [comp_cost] = cost_function(avg_total_intensity, T, m, h)
    start_path = getsampleusingtime(avg_total_intensity, 0, T-h);
    end_path = getsampleusingtime(avg_total_intensity, T-h, T);

    if length(start_path.Time) == 1
        start_lambda_integral = (T-h) * start_path.Data;
    else
        start_lambda_integral = trapz(start_path.Time, start_path.Data);
    end

    if length(end_path.Time) == 1
        end_lambda_integral = h * end_path.Data;
    else
        end_lambda_integral = trapz(end_path.Time, end_path.Data);
    end

    comp_cost = start_lambda_integral + end_lambda_integral * m;
end

