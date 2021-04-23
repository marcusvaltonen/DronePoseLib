function [f_err, F_err, r_err] = compare_to_gt(out, f, F, r)

N = length(out);

tmp1 = zeros(N, 1);
tmp2 = zeros(N, 1);

if nargin < 4
    tmp3 = zeros(N, 1);
end

% Normalize GT
F = F / F(3, 3);

for k = 1:N
    f_est = out(k).f;
    tmp1(k) = abs(f_est - f) / f;
    F_est = out(k).F;
    tmp2(k) = norm(F_est / F_est(3, 3) - F) / norm(F);
    if nargin > 3
        r_est = out(k).r;
        tmp3(k) = abs(r_est - r) / abs(r);
    end
end

% Compute errors
f_err = min(tmp1);
F_err = min(tmp2);
if nargout > 2
    r_err = min(tmp3);
end