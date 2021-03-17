% Script to compile and link C++ source file to MATLAB compatible MEX binary.
solver_types = {'get_fEf', 'get_frEfr'};
path_to_eigen = '/usr/include/eigen3';
cpp_flags = '-O';
outputdir = 'solvers';

if ~exist(outputdir, 'dir')
    mkdir(outputdir)
end

for j = 1:length(solver_types)
    switch solver_types{j}
        case 'get_fEf'
            mex(['-I', path_to_eigen], ...
                '-I../includes/DronePoseLib', ...
                '-I../src/helpers', ...
                'src/get_valtonenornhag_arxiv_2021_fEf_mex.cpp', ...
                '../src/solvers/valtonenornhag_arxiv_2021/fEf/get_fEf.cpp', ...
                '../src/helpers/normalize2dpts.cpp', ...
                '../src/helpers/quartic.cpp', ...
                cpp_flags, '-outdir', outputdir)
        case 'get_frEfr'
            mex(['-I', path_to_eigen], ...
                '-I../includes/DronePoseLib', ...
                '-I../src/helpers', ...                
                '-I../src/solvers/valtonenornhag_arxiv_2021/frEfr', ...
                'src/get_valtonenornhag_arxiv_2021_frEfr_mex.cpp', ...
                '../src/helpers/normalize2dpts.cpp', ...
                '../src/helpers/radial.cpp', ...
                '../src/solvers/valtonenornhag_arxiv_2021/frEfr/coeffs_frEfr.cpp', ...
                '../src/solvers/valtonenornhag_arxiv_2021/frEfr/get_frEfr.cpp', ...
                '../src/solvers/valtonenornhag_arxiv_2021/frEfr/solver_frEfr.cpp', ...
                cpp_flags, '-outdir', outputdir)
        otherwise
            error('Not a valid method.')
    end
end
