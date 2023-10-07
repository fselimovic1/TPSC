function [ results, measurements, data ] = runsse(casename, vrs, sesettings)
% ----------------- Load Power System Data and Measurements ---------------
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');
% -------------------------------------------------------------------------

%----------------- Extract Useful Informations (Power System) -------------
powsys = preprocess_ps(data, 'se');
%--------------------------------------------------------------------------

% --------------------- Calculate Y matrix --------------------------------
powsys = admittance_matrix(powsys);
% -------------------------------------------------------------------------

%----------------- Extract Useful Informations (Measurements) -------------
meas = preprocess_meas(powsys, measurements);
%--------------------------------------------------------------------------

% ---------------- Solve the state estimation problem ---------------------
tic
if strcmp('complex', sesettings.domain)
    if strcmp('cgn_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, sesettings);
    elseif strcmp('ecgn_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_ecgn_sse(powsys, meas, sesettings);
    elseif strcmp('cec_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_cec_sse(powsys, meas, sesettings);
    elseif strcmp('pgne', sesettings.method)
        sesettings.initialStage = 1;
        [ Vc, iter, converged, info ] = run_pgne_rtse(powsys, meas, sesettings, ones(powsys.num.bus, 1));
    elseif strcmp('cls_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_cls_sse(powsys, meas);
    elseif strcmp('lcec_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_lcec_sse(powsys, meas);
    else
        error("Specified method does not exist.")
    end
else
    if strcmp('wls_rect_sse', sesettings.method)
        [ Vc, iter, converged, info ] = run_wls_rect_sse(powsys, meas);
        Vc = Vc(1:2:2*powsys.num.bus - 1) + 1i * Vc(2:2:2*powsys.num.bus);
    else
    end
end
algtime = toc;
% -------------------------------------------------------------------------

% -------------------------- Show results ---------------------------------
results.voltage = Vc(powsys.bus.busnew);
results.converged = converged;
results.iter = iter;
results.algtime = algtime;
results.info = info;
results.sys = casename;

results_sse(results, powsys, meas, sesettings)
% -------------------------------------------------------------------------
end

