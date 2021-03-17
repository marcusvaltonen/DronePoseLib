addpath('helpers')
addpath(genpath('solvers'))

nbr_iter = 1e4;
error_F = zeros(nbr_iter,1);
error_f = zeros(nbr_iter,1);

% This is a 3-point method
N = 3;

for j = 1:nbr_iter
    % Generate random problem instance
    [R1, R2, Ry, f, F, x1, x2, R, t, x1u, x2u] = generate_points_realistic(N, 0, 0);
    
    % Solve the problem
    out = get_valtonenornhag_arxiv_2021_fEf_mex(x1(1:2,:), x2(1:2,:), R1, R2);
    
    % Compare to ground truth
    [error_f(j), error_F(j)] = compare_to_gt(out, f, F);
end

%% Plot histogram
figure(1)
histogram(log10(error_f), 50)
title('Focal length error')

figure(2)
histogram(log10(error_F), 50)
title('Fundamental matrix error')