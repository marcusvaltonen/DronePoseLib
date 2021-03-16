function out = get_fEf_no_rot_wrapper(x1, x2, R1, R2)
% The MEX-file does not compute the fundamental matrix.
% TODO: Include these computations in the C++ file.

out = get_valtonenornhag_arxiv_2021_fEf_mex(x1, x2, R1, R2);

for j = 1:length(out)
    f = out(j).f;
    Kinv = diag([1 / f, 1 / f, 1]);
    R = R2 * R1';
    t = out(j).t;
    out(j).F = Kinv * skew(t) * R * Kinv;
end