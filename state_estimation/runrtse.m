function [ results, measurements, data ] = runrtse(casename, vrs, rtsesettings)
% --------------------- Load PS Data and Measurements ---------------------
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');
% -------------------------------------------------------------------------

% ------------------------ Settings Check ---------------------------------
if ~strcmp(measurements.mode, "tracking")
	error('Measurements are not generated in tracking mode. Real time state estimation reqiures measurement for a time interval.');
end
if measurements.genFreq < rtsesettings.fc
	error('Available measurements does not support the entered estimation frequency.');
end
dI = measurements.genFreq/rtsesettings.fc;
if dI ~= fix(dI)
	error('Estimation/tracking frequency must be an integer divisor of the base measurements recieving frequency.');
end
estamps = ceil(measurements.tstamps/dI);
% -------------------------------------------------------------------------

%----------------- Extract Useful Informations (Power System) -------------
powsys = preprocess_ps(data, 'se');
%--------------------------------------------------------------------------

% --------------------- Calculate Y matrix --------------------------------
powsys = admittance_matrix(powsys);
% -------------------------------------------------------------------------

% -------------------- Allocate memeory for the results -------------------
results.Vc = complex(zeros(powsys.num.bus, estamps));
% -------------------------------------------------------------------------

% ---------------------- Prepare real time plot ---------------------------
if rtsesettings.realtimeplot
    rtpEst = animatedline('Color','r');
    rtpTrue = animatedline('Color','k', 'LineWidth', 1.5);
%     axis([0 round(measurements.tstamps/measurements.genFreq) 0.9 1.1])
    str = ['Bus ', num2str(rtsesettings.rtpbus), ' voltage magnitude'];
    title(str)
    xlabel("Time [s]")
    ylabel("Voltage magnitude [p.u.]")
    legend("Estimated", "True")
end
% -------------------------------------------------------------------------
for i = 1:dI:measurements.tstamps
    % -------------------- Extract Present Measurements -------------------
    cMeasurements.scada = measurements.scada(measurements.scada(:, 1) == i - 1, :);
	cMeasurements.synpmu = measurements.synpmu(measurements.synpmu(:, 1) == i - 1, :);
    cMeasurements.trueVoltage = measurements.trueVoltage(:, i);
    cMeasurements.fpmu = measurements.fpmu(measurements.fpmu(:, 1) == i - 1, :);
    cMeasurements.mode = 'static';
    meas = preprocess_meas(powsys, cMeasurements);
    % ---------------------------------------------------------------------
    
    % ------------------- Calculate state variables -----------------------
    if strcmp(rtsesettings.domain, "complex")
        if strcmp(rtsesettings.method, "quasidyn") 
            % --------------------- Initial state variables ---------------
            if rtsesettings.initialStage
                if rtsesettings.flatStart     
                    Vc = ones(2  * powsys.num.bus, 1);
                else
                    Vc = [ powsys.bus.Vmi .* exp(1i * powsys.bus.Vai);
                           powsys.bus.Vmi .* exp(-1i * powsys.bus.Vai); 
                           ];
                end
            end
            % -------------------------------------------------------------
            [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, rtsesettings, Vc);
        else
        end
    else
        if strcmp(rtsesettings.method, 'quasidyn')
                [ X, iter, converged, info  ] = run_wls_rect_sse(powsys,...
                meas);
        elseif strcmp(rtsesettings.method, 'fEKFrect')
            % --------------------- Initial state variables ---------------
            if rtsesettings.initialStage
                if rtsesettings.flatStart     
                    X_ = [ reshape([ones(1, powsys.num.bus, 1); ...
                           zeros(1 ,powsys.num.bus) ], [], 1); 0 ];
                else
                    X_ = [ powsys.bus.Vmi .* cos(powsys.bus.Vai);
                           powsys.bus.Vmi .* sin(powsys.bus.Vai); 
                           0
                           ];
                end
            end
            % -------------------------------------------------------------
            [ X, X_, converged, info ] = run_fEKFrect_rtse(powsys,...
                meas, rtsesettings, X_);                
        end
        Vc = X(1:2:2 * powsys.num.bus - 1) + 1i .* X(2:2:2 * powsys.num.bus);
    end
    % ---------------------------------------------------------------------
    if ~converged
        error("State estimation algorithm did not converged. Real-time state estimation calculations stopped.");
    end
    % ---------------------- Real time plot -------------------------------
    if rtsesettings.realtimeplot
        addpoints(rtpEst, (i-1)/measurements.genFreq, abs(Vc(rtsesettings.rtpbus)));
        addpoints(rtpTrue, (i-1)/measurements.genFreq, abs(cMeasurements.trueVoltage(rtsesettings.rtpbus)));
        drawnow;
    end
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    pause(0.1)
    rtsesettings.initialStage = false;
end
end

