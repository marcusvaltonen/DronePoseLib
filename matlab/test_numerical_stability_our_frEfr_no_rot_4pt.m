addpath('helpers')
addpath(genpath('solvers'))

nbr_iter = 1e3;
error_f = zeros(nbr_iter,1);
error_F = zeros(nbr_iter,1);
error_r = zeros(nbr_iter,1);

% This is a 4-point method
N = 4;

for j = 1:nbr_iter
    % Generate random problem instance
    r = -rand * (1 - 1e-6);
    [R1, R2, ~, f, F, x1, x2] = generate_points(N, 0, r);

    % Solve the problem
    out = get_frEfr_no_rot_wrapper(x1(1:2,:), x2(1:2,:), R1, R2);

    % Compare to ground truth
    tmp1 = [];
    tmp2 = [];
    tmp3 = [];

    % Normalize GT
    F = F / F(3, 3);
    
    for k = 1:length(out)
        f_est = out(k).f;
        tmp1(k) = abs(f_est - f) / f;
        r_est = out(k).r;
        tmp2(k) = abs(r_est - r) / abs(r);
        F_est = out(k).F / out(k).F(3, 3);
        tmp4(k) = norm(F_est - F) / norm(F);
    end

    
    % Compute errors
    error_f(j) = min(tmp1);
    error_r(j) = min(tmp2);
    error_F(j) = min(tmp4);
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