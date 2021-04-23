addpath('helpers')
addpath(genpath('solvers'))

nbr_iter = 1e4;
error_f = zeros(nbr_iter,1);
error_F = zeros(nbr_iter,1);
error_r = zeros(nbr_iter,1);

% This is a 4-point method
N = 4;

for j = 1:nbr_iter
    if mod(j, 500) == 0
        fprintf(1, 'Iter: %d\n', j)
    end
    % Generate random problem instance
    r = -rand * 1e-6;
    [R1, R2, ~, f, F, x1, x2, R, t, x1u, x2u] = generate_points_realistic(N, 0, r);

    % Solve problem
    out = get_valtonenornhag_arxiv_2021_frEfr_mex(x1(1:2,:), x2(1:2,:), R1, R2);

    % Compare to ground truth
    [error_f(j), error_F(j), error_r(j)] = compare_to_gt(out, f, F, r);
end

%% Plot histogram
figure(1)
histogram(log10(error_f), 50)
title('Focal length error')

figure(2)
histogram(log10(error_F), 50)
title('Fundamental matrix error')

figure(3)
histogram(log10(error_r), 50)
title('Radial distortion parameter error')