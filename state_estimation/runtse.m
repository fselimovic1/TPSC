function [ results, measurements, data ] = runtse(casename, vrs, tsesettings)
% --------------------- Load PS Data and Measurements ---------------------
data = loadcase(casename, vrs);
measurements = loadcase(casename, vrs, 'M');
% -------------------------------------------------------------------------

% ------------------------ Settings Check ---------------------------------
if ~strcmp(measurements.mode, "tracking")
	error('Measurements are not generated in tracking mode. Real time state estimation reqiures measurement for a time interval.');
end
if measurements.genFreq < tsesettings.fc
	error('Available measurements do not support the entered estimation frequency.');
end
dI = measurements.genFreq/tsesettings.fc;
if dI ~= fix(dI)
	error('Estimation/tracking frequency must be an integer divisor of the base measurements recieving frequency.');
end
tstamps = ceil(measurements.tstamps/dI);
% -------------------------------------------------------------------------

%----------------- Extract Useful Informations (Power System) -------------
powsys = preprocess_ps(data, 'se');
%--------------------------------------------------------------------------

% --------------------- Calculate Y matrix --------------------------------
powsys = admittance_matrix(powsys);
% -------------------------------------------------------------------------

% -------------------- Allocate memeory for the results -------------------
global X
results.Vc = complex(zeros(powsys.num.bus, tstamps));
rmse.vm = zeros(tstamps, 1);
rmse.va = zeros(tstamps, 1);
if strcmp(tsesettings.method, 'fEKFrect')
    results.f = zeros(tstamps, 1);
end
% -------------------------------------------------------------------------


% ------------------ OPTIONAL - Measure execution time --------------------
if tsesettings.measureTime
    times  = zeros(tstamps - 1, 1); 
end
% -------------------------------------------------------------------------

% ---------------------- Prepare real time plot ---------------------------
if tsesettings.realtimeplot
    % Define the figure size in inches (width and height)
    width = 3.5;   % Width of the figure in inches
    height = 2.8;  % Height of the figure in inches
    % Create a new figure with the specified size
    figure('Position', [100, 100, width * 100, height * 100]); % Position is in pixels

    rtpEst = animatedline('color', 'r');
    rtpTrue = animatedline('Color', 'k', 'LineWidth', 1.5);
%     axis([0 round(measurements.tstamps/measurements.genFreq) 0.9 1.1])
    if ~tsesettings.plotForPaper
        str = ['Bus ', num2str(tsesettings.rtpbus), ' voltage magnitude'];
        title(str)
        xlabel("Time [s]")
        ylabel("Voltage magnitude [p.u.]")
        legend("Estimated", "True")
    end
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
    tic
    % ---------------------- COMPLEX DOMAIN METHODS -----------------------
    % ---------------------------------------------------------------------
    if strcmp(tsesettings.domain, "complex")
        if strcmp(tsesettings.method, "quasidyn") 
            % --------------------- Initial state variables ---------------
            if tsesettings.initialStage
                if tsesettings.flatStart     
                    Vc = ones(2  * powsys.num.bus, 1);
                else
                    Vc = [ powsys.bus.Vmi .* exp(1i * powsys.bus.Vai);
                           powsys.bus.Vmi .* exp(-1i * powsys.bus.Vai); 
                           ];
                end
            end
            % -------------------------------------------------------------
            [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, tsesettings, Vc);
        elseif strcmp(tsesettings.method, "pgne")
            % --------------------- Initial state variables ---------------
            if tsesettings.initialStage
                if tsesettings.flatStart     
                    Vc = ones(powsys.num.bus, 1);
                else
                    Vc = powsys.bus.Vmi .* exp(1i * powsys.bus.Vai);
                end
            end
            [ Vc, iter, converged, info ] = run_pgne_tse(powsys, meas, tsesettings, Vc);
            % -------------------------------------------------------------
        elseif strcmp(tsesettings.method, "ckf")
            %--------------------- Initial state variables ----------------
            if tsesettings.initialStage
                if tsesettings.flatStart     
                    x_ = ones(powsys.num.bus, 1);
                else
                    x_ = powsys.bus.Vmi .* (cos(powsys.bus.Vai)...
                            +  1i .* sin(powsys.bus.Vai)); 
                end
            end
            % -------------------------------------------------------------
            if tsesettings.initialStage
                X = zeros(powsys.num.bus, tstamps);
            end
            tsesettings.tstep = ceil(i/dI);
            [ Vc, x_, converged, info ] = run_ckf_dse(powsys,...
                meas, tsesettings, x_);
            X(:, ceil(i/dI)) = Vc;
        end
    else
    % ------------------------ REAL DOMAIN METHODS ------------------------
    % ---------------------------------------------------------------------
        if strcmp(tsesettings.method, 'quasidyn')
                [ x, iter, converged, info  ] = run_wls_rect_sse(powsys,...
                meas, tsesettings);
        elseif strcmp(tsesettings.method, 'fEKFrect')
            %--------------------- Initial state variables ----------------
            if tsesettings.initialStage
                if tsesettings.flatStart     
                    x_ = [ reshape([ones(1, powsys.num.bus, 1); ...
                           zeros(1 ,powsys.num.bus) ], [], 1); 0 ];
                else
                    x_ = [ powsys.bus.Vmi .* cos(powsys.bus.Vai);
                           powsys.bus.Vmi .* sin(powsys.bus.Vai); 
                           0
                           ];
                end
            end
            % -------------------------------------------------------------
            [ x, x_, converged, info ] = run_fEKFrect_dse(powsys,...
                meas, tsesettings, x_);       
            results.f(i) = 50 + x(2 * powsys.num.bus + 1);
        elseif strcmp(tsesettings.method, "qDKF")
            %--------------------- Initial state variables ----------------
            if tsesettings.initialStage
                if tsesettings.flatStart     
                    x_ = [ reshape([ones(1, powsys.num.bus, 1); ...
                           zeros(1 ,powsys.num.bus) ], [], 1) ];
                else
                    x_ = [ powsys.bus.Vmi .* cos(powsys.bus.Vai);
                           powsys.bus.Vmi .* sin(powsys.bus.Vai); 
                           ];
                end
            end
            % -------------------------------------------------------------
            if tsesettings.initialStage
                X = zeros( 2 * powsys.num.bus, tstamps);
            end
            tsesettings.tstep = ceil(i/dI);
            [ x, x_, converged, info ] = run_qDKF_dse(powsys,...
                meas, tsesettings, x_);
            X(:, ceil(i/dI)) = x;
        end
        Vc = x(1:2:2 * powsys.num.bus - 1) + 1i .* x(2:2:2 * powsys.num.bus);
    end
    % ------------------ OPTIONAL - Measure execution time ----------------
    if tsesettings.measureTime && i ~= 1
        times(ceil(i/dI)) = toc;
    end
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    if ~converged
        error("State estimation algorithm did not converged. Real-time state estimation calculations stopped.");
    end
    % ---------------------- Real time plot -------------------------------
    Vc = full(Vc);
    if tsesettings.realtimeplot
        addpoints(rtpEst, (i-1)/measurements.genFreq, abs(Vc(tsesettings.rtpbus)));
        addpoints(rtpTrue, (i-1)/measurements.genFreq, abs(cMeasurements.trueVoltage(tsesettings.rtpbus)));
        drawnow;
    end
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    % ----------------------- Compute RMSE  -------------------------------
    rmse.vm(ceil(i/dI)) = sqrt(sum(abs(abs(measurements.trueVoltage(:, i))... 
                        - abs(Vc(1:powsys.num.bus)))));
    rmse.va(ceil(i/dI)) = sqrt(sum(abs(angle(measurements.trueVoltage(:, i))...
                        - angle(Vc(1:powsys.num.bus)))));
    % ---------------------------------------------------------------------
    results.Vc(:, i) = Vc;
    pause(tsesettings.plotpause)
    tsesettings.initialStage = false;
end
if tsesettings.measureTime
    fprintf("%s requires %.2f [ms] in average.\n", tsesettings.method, mean(times) * 1000);
end
results.rmse = rmse;
end

