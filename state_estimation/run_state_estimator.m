function [results] = run_state_estimator(sesettings, data, measurements)
if strcmp('complex', sesettings.domain)
    results = run_cgn_sse(sesettings, data, measurements);
elseif strcmp('real', sesettings.domain)
    if strcmp('sgn_sse', sesettings.method)
        results = run_gn_sse(sesettings, data, measurements);
    elseif strcmp('wls_tse', sesettings.method)
        results = run_wls_tse_pmu(sesettings, data, measurements);
    end
end
end

