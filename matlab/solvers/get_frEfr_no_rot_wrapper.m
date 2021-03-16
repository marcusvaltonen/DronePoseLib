function out = get_frEfr_no_rot_wrapper(x1, x2, R1, R2)
% The MEX-file does not compute the fundamental matrix, nor normalization.
% TODO: Include these computations in the C++ file.

assert(all(size(x1) == [2, 4]))
assert(all(size(x2) == [2, 4]))
assert(all(size(R1) == [3, 3]))
assert(all(size(R2) == [3, 3]))

% Normalize
x1 = [x1; ones(1, 4)];
x2 = [x2; ones(1, 4)];
scale = max(normalise2dpts(x1), normalise2dpts(x2));
S = [scale 0 0; 0 scale 0; 0 0 1];
x1 = pflat(S*x1);
x2 = pflat(S*x2);

% Run C++ solver
out = get_valtonenornhag_arxiv_2021_frEfr_mex(x1(1:2, :), x2(1:2, :), R1, R2);

for j = 1:length(out)
    out(j).f = out(j).f / scale;
    out(j).r = out(j).r * scale^2;
    f = out(j).f;
    Kinv = diag([1 / f, 1 / f, 1]);
    R = R2 * R1';
    t = out(j).t;
    out(j).F = Kinv * skew(t) * R * Kinv;
end