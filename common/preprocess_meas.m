function meas = preprocess_meas(data, measurements)
% ------------------------- SCADA Measurements ----------------------------
if ~isempty(measurements.scada) && any(measurements.scada(:, 1) ~= measurements.scada(1, 1))
    error('When runnning static state estimation, all measurements must have the same time stamp.')
end
meas.scada.type = measurements.scada(:, 3);
meas.scada.onbus = measurements.scada(:, 2);
meas.scada.loc = measurements.scada(:, 4);
meas.scada.m = measurements.scada(:, 5);
meas.scada.true = measurements.scada(:, 6);
meas.scada.pbranchO = find(meas.scada.type == 1 & meas.scada.loc < 0);
meas.scada.pbranch = find(meas.scada.type == 1 & meas.scada.loc > 0);
meas.scada.qbranchO = find(meas.scada.type == 2 & meas.scada.loc < 0);
meas.scada.qbranch = find(meas.scada.type == 2 & meas.scada.loc > 0);
meas.scada.pinj = find(meas.scada.type == 3);
meas.scada.qinj = find(meas.scada.type == 4);
meas.scada.ibrMO = find(meas.scada.type == 5);
meas.scada.ibrM = find(meas.scada.type == 5 & meas.scada.loc < 0);
meas.scada.vm = find(meas.scada.type == 6);
% meas.scada.sd = data.scada(:, 4);
% -------------------------------------------------------------------------

% -------------------------- PMU Measurements -----------------------------
if ~isempty(measurements.scada) && any(measurements.synpmu(:, 1) ~= measurements.synpmu(1, 1))
    error('When runnning static state estimation, all measurements must have the same time stamp.')
end
meas.pmu.type = measurements.synpmu(:, 3);
meas.pmu.onbus = measurements.synpmu(:, 2);
meas.pmu.loc = measurements.synpmu(:, 4);
meas.pmu.m = measurements.synpmu(:, 5);
meas.pmu.a = measurements.synpmu(:, 6);
meas.pmu.mtrue = measurements.synpmu(:, 7);
meas.pmu.atrue = measurements.synpmu(:, 8);
% meas.pmu.msd = data.pmu(meas.pmu.onbus, 3); 
% meas.pmu.asd = data.pmu(meas.pmu.onbus, 4);
meas.pmu.inj = find(meas.pmu.type == 2);
meas.pmu.ibranch = find(meas.pmu.type == 1 & meas.pmu.loc > 0);
meas.pmu.ibranchO = find(meas.pmu.type == 1 & meas.pmu.loc < 0);
meas.pmu.vnode = find(meas.pmu.type == 3);
% -------------------------------------------------------------------------

% ------------------------ Number of measurements -------------------------
meas.num.scada = size(measurements.scada, 1);
meas.num.pmu = size(measurements.synpmu, 1);

% ----------------------------- SCADA (-s) --------------------------------
meas.num.sPbr0 = numel(meas.scada.pbranchO);
meas.num.sPbr = numel(meas.scada.pbranch);
meas.num.sQbr0 = numel(meas.scada.qbranchO);
meas.num.sQbr = numel(meas.scada.qbranch);
meas.num.sPinj = numel(meas.scada.pinj);
meas.num.sQinj = numel(meas.scada.qinj);
meas.num.sIbrM0 = numel(meas.scada.ibrMO);
meas.num.sIbrM = numel(meas.scada.ibrM);
meas.num.sVm = numel(meas.scada.vm);
% -------------------------------------------------------------------------

% ------------------------------ PMU (-p) ---------------------------------
meas.num.pInj = numel(meas.pmu.inj);
meas.num.pIij = numel(meas.pmu.ibranch);
meas.num.pIijO = numel(meas.pmu.ibranchO);
meas.num.pV = numel(meas.pmu.vnode);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% --------------------------- True values ---------------------------------
meas.Vtrue = measurements.trueVoltage;
% -------------------------------------------------------------------------
end

