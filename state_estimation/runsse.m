function [ results, measurements, data ] = runsse(casename, vrs, sesettings)
% ----------------- Load Power System Data and Measurements ---------------
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');
% -------------------------------------------------------------------------

%--------------------- Extract Useful Informations ------------------------
powsys = preprocess(data, 'se');
%--------------------------------------------------------------------------


if strcmp('complex', sesettings.domain)
    results = run_cgn_sse(sesettings, data, measurements);
elseif strcmp('real', sesettings.domain)
    if strcmp('sgn_sse', sesettings.method)
        results = run_gn_sse(sesettings, data, measurements);
    elseif strcmp('wls_tse', sesettings.method)
        results = run_wls_tse_pmu(sesettings, data, measurements);
    else
        return
    end
else
    return
end
results_sse(sesettings, measurements, results)
end

