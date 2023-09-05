function meas = preprocess_meas(powsys, measurements)
if ~strcmp(measurements.mode, "static")
    error('Measurements are not generated in static mode. Static state estimation reqiures measurement for a distinct time stamp.')
end
% ------------------------- SCADA Measurements ----------------------------
meas.scada.type = measurements.scada(:, 3);
meas.scada.onbus = measurements.scada(:, 2);
meas.scada.loc = measurements.scada(:, 4);
meas.scada.m = measurements.scada(:, 5);
meas.scada.true = measurements.scada(:, 6);
meas.scada.pijO = find(meas.scada.type == 1 & meas.scada.loc < 0);
meas.scada.pij = find(meas.scada.type == 1 & meas.scada.loc > 0);
meas.scada.qijO = find(meas.scada.type == 2 & meas.scada.loc < 0);
meas.scada.qij = find(meas.scada.type == 2 & meas.scada.loc > 0);
meas.scada.pinj = find(meas.scada.type == 3);
meas.scada.qinj = find(meas.scada.type == 4);
meas.scada.IijmO = find(meas.scada.type == 5 & meas.scada.loc < 0);
meas.scada.Iijm = find(meas.scada.type == 5 & meas.scada.loc > 0);
meas.scada.vm = find(meas.scada.type == 6);
% -------------------------------------------------------------------------

% -------------------------- PMU Measurements -----------------------------
meas.pmu.type = measurements.synpmu(:, 3);
meas.pmu.onbus = measurements.synpmu(:, 2);
meas.pmu.loc = measurements.synpmu(:, 4);
meas.pmu.m = measurements.synpmu(:, 5);
meas.pmu.a = measurements.synpmu(:, 6);
meas.pmu.mtrue = measurements.synpmu(:, 7);
meas.pmu.atrue = measurements.synpmu(:, 8);
meas.pmu.Iinj = find(meas.pmu.type == 2);
meas.pmu.Iij = find(meas.pmu.type == 1 & meas.pmu.loc > 0);
meas.pmu.IijO = find(meas.pmu.type == 1 & meas.pmu.loc < 0);
meas.pmu.v = find(meas.pmu.type == 3);

% --------------------- Frequency measurements ----------------------------
if isfield(measurements, 'fpmu')
    meas.fpmu.onbus = measurements.fpmu(:, 2);
    meas.fpmu.m = measurements.fpmu(:, 3);
    meas.num.f = numel(meas.fpmu.onbus);
end
% -------------------------------------------------------------------------
meas.pmu.msd = ((powsys.pmu.msd(meas.pmu.onbus) .* meas.pmu.m)/100) / 3;
meas.pmu.asd = (powsys.pmu.asd(meas.pmu.onbus) * pi / 180)/3;
meas.fpmu.fsd = (powsys.pmu.fsd(meas.fpmu.onbus)) / 3;
% -------------------------------------------------------------------------

% ------------------------ Number of measurements -------------------------
meas.num.scada = size(measurements.scada, 1);
meas.num.pmu = size(measurements.synpmu, 1);

% ----------------------------- SCADA (-s) --------------------------------
meas.num.sPijO = numel(meas.scada.pijO);
meas.num.sPij = numel(meas.scada.pij);
meas.num.sQijO = numel(meas.scada.qijO);
meas.num.sQij = numel(meas.scada.qij);
meas.num.sPinj = numel(meas.scada.pinj);
meas.num.sQinj = numel(meas.scada.qinj);
meas.num.sIijmO = numel(meas.scada.IijmO);
meas.num.sIijm = numel(meas.scada.Iijm);
meas.num.sVm = numel(meas.scada.vm);
% -------------------------------------------------------------------------

% ------------------------------ PMU (-p) ---------------------------------
meas.num.pIinj = numel(meas.pmu.Iinj);
meas.num.pIij = numel(meas.pmu.Iij);
meas.num.pIijO = numel(meas.pmu.IijO);
meas.num.pV = numel(meas.pmu.v);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% --------------------------- True values ---------------------------------
meas.Vtrue = measurements.trueVoltage;
% -------------------------------------------------------------------------
end

