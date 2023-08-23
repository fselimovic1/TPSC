function [ results, measurements, data ] = run_state_estimator(casename, vrs, sesettings)
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');
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

