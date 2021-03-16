function [f_err, F_Err] = compare_to_gt(out, f, F)

N = length(out);

tmp1 = zeros(N,1);
tmp2 = zeros(N,1);

% Normalize GT
F = F / F(3, 3);

for k = 1:N
    f_est = out(k).f;
    tmp1(k) = abs(f_est - f) / f;
    F_est = out(k).F;
    tmp2(k) = norm(F_est / F_est(3, 3) - F) / norm(F);
end

% Compute errors
f_err = min(tmp1);
F_Err = min(tmp2);