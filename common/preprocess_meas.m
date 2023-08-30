function meas = preprocess_meas(data, measurements)
% ------------------------- SCADA Measurements ----------------------------
if any(measurements.scada(:, 1) ~= measurements.scada(1, 1))
    error('When runnning static state estimation, all measurements must have the same time stamp.')
end
meas.scada.type = measurements.scada(:, 3);
meas.scada.onbus = measurements.scada(:, 2);
meas.scada.loc = measurements.scada(:, 4);
meas.scada.m = measurements.scada(:, 5);
meas.scada.true = measurements.scada(:, 6);
meas.scada.sd = data.scada(:, 4);
% -------------------------------------------------------------------------

% -------------------------- PMU Measurements -----------------------------
if any(measurements.synpmu(:, 1) ~= measurements.synpmu(1, 1))
    error('When runnning static state estimation, all measurements must have the same time stamp.')
end
meas.pmu.type = measurements.pmu(:, 3);
meas.pmu.onbus = measurements.pmu(:, 2);
meas.pmu.loc = measurements.pmu(:, 4);
meas.pmu.m = measurements.scada(:, 5);
meas.pmu.a = measurements.scada(:, 6);
meas.pmu.mtrue = measurements.scada(:, 7);
meas.pmu.atrue = measurements.scada(:, 8);
meas.pmu.msd = data.pmu(meas.pmu.onbus, 3); 
meas.pmu.asd = data.pmu(meas.pmu.onbus, 4);
meas.pmu.ibranch = find(meas.pmu.type == 1);
meas.pmu.vnode = find(meas.pmu.type == 3);
% -------------------------------------------------------------------------

% ------------------------ Number of measurements -------------------------
meas.num.scada = size(measurements.scada, 1);
meas.num.pmu = size(measurements.synpmu, 1);
meas.num.pIij = numel(meas.pmu.ibranch);
meas.num.pV = numel(meas.pmu.vnode);
% -------------------------------------------------------------------------
end

