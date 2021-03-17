addpath('helpers')
addpath(genpath('solvers'))

nbr_iter = 100;
noise_levels = 1:5;
nbr_noise_levels = length(noise_levels);

error_F_3pt = nan(nbr_noise_levels, nbr_iter);
error_f_3pt = nan(nbr_noise_levels, nbr_iter);
error_F_4pt = nan(nbr_noise_levels, nbr_iter);
error_f_4pt = nan(nbr_noise_levels, nbr_iter);

N = 4;

for i = 1:nbr_noise_levels
    disp(['Noise level: ', num2str(i)])
    noise = noise_levels(i);
    for j = 1:nbr_iter
        % Add noise on R1 and R2.
        R1noise = eul2rotm(deg2rad(0.1 * noise * randn(1, 3)), 'XYZ');
        R2noise = eul2rotm(deg2rad(0.1 * noise * randn(1, 3)), 'XYZ');

        % Generate points on ground plane
        [R1_orig, R2_orig, Ry, f, F, x1, x2, R, t, x1u, x2u] = generate_points_realistic(N, 0, 0);

        R1 = R1_orig * R1noise;
        R2 = R2_orig * R2noise;
        
        % Add noise to rectified image points
        p1 = [x1(1:2, :) + f / 1080 * randn(2, N); ones(1, N)];
        p2 = [x2(1:2, :) + f / 1080 * randn(2, N); ones(1, N)];
        
        % Valtonen Ornhag etal. (Arxiv, 2021) - fEf (3 pt)
        try
            out = get_valtonenornhag_arxiv_2021_fEf_mex(p1(1:2, 1:3), p2(1:2, 1:3), R1, R2);
            [error_f_3pt(i,j), error_F_3pt(i,j)] = compare_to_gt(out, f, F);
        catch
            warning('3pt failed')
        end

        % Add noise to unrectified image points
        p1 = [x1u(1:2, :) + f / 1080 * randn(2, N); ones(1, N)];
        p2 = [x2u(1:2, :) + f / 1080 * randn(2, N); ones(1, N)];
        
        % Valtonen Ornhag etal. (Arxiv, 2021) - frEfr (4 pt)
        try
            out = get_valtonenornhag_arxiv_2021_frEfr_mex(p1(1:2, 1:4), p2(1:2, 1:4), R1, R2);
            [error_f_4pt(i,j), error_F_4pt(i,j)] = compare_to_gt(out, f, F);
        catch
            warning('4pt failed')
        end        
         
    end
end
%% Draw boxplots
figure(1)
subplot(121)
boxplot(log10(error_F_3pt'))
xlabel('3pt')
ylim([-6, 0])
subplot(122)
boxplot(log10(error_F_4pt'))
xlabel('4pt')
ylim([-6, 0])
sgtitle('Fundamental matrix error')

figure(2)
subplot(121)
boxplot(log10(error_f_3pt'))
xlabel('3pt')
ylim([-5, 1])
subplot(122)
boxplot(log10(error_f_4pt'))
xlabel('4pt')
ylim([-5, 1])
sgtitle('Focal length error')