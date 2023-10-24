function meas = preprocess_meas(powsys, measurements)
if ~strcmp(measurements.mode, "static")
    error(['Measurements are not generated in static mode.',...
        'Static state estimation reqiures measurement for a distinct time stamp.'])
end
% ------------------------- SCADA Measurements ----------------------------
meas.scada.type = measurements.scada(:, 3);
meas.scada.onbus = measurements.scada(:, 2);
meas.scada.loc = measurements.scada(:, 4);
meas.scada.m = measurements.scada(:, 5);
meas.scada.true = measurements.scada(:, 6);
meas.scada.pji = find(meas.scada.type == 1 & meas.scada.loc < 0);
meas.scada.pij = find(meas.scada.type == 1 & meas.scada.loc > 0);
meas.scada.qji = find(meas.scada.type == 2 & meas.scada.loc < 0);
meas.scada.qij = find(meas.scada.type == 2 & meas.scada.loc > 0);
meas.scada.pi = find(meas.scada.type == 3);
meas.scada.qi = find(meas.scada.type == 4);
meas.scada.Ijim = find(meas.scada.type == 5 & meas.scada.loc < 0);
meas.scada.Iijm = find(meas.scada.type == 5 & meas.scada.loc > 0);
meas.scada.vm = find(meas.scada.type == 6);
% ---------------- Complex pairs (needed for some estimators) -------------
[ meas.scada.sji.p, meas.scada.sji.q ] = myintersect(meas.scada.loc, meas.scada.pji, meas.scada.qji);
[ meas.scada.sij.p, meas.scada.sij.q ] = myintersect(meas.scada.loc, meas.scada.pij, meas.scada.qij);
[ meas.scada.si.p, meas.scada.si.q ] = myintersect(meas.scada.loc, meas.scada.pi, meas.scada.qi);
% -------------------------------------------------------------------------

% ----------------------- Real powers indices only ------------------------
meas.scada.oddPji = setdiff(meas.scada.pji, meas.scada.sji.p);
meas.scada.oddQji = setdiff(meas.scada.qji, meas.scada.sji.q);
meas.scada.oddPij = setdiff(meas.scada.pij, meas.scada.sij.p);
meas.scada.oddQij = setdiff(meas.scada.qij, meas.scada.sij.q);
meas.scada.oddPi = setdiff(meas.scada.pi, meas.scada.si.p);
meas.scada.oddQi = setdiff(meas.scada.qi, meas.scada.si.q);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------- PMU Measurements -----------------------------
meas.pmu.type = measurements.synpmu(:, 3);
meas.pmu.onbus = measurements.synpmu(:, 2);
meas.pmu.loc = measurements.synpmu(:, 4);
meas.pmu.m = measurements.synpmu(:, 5);
meas.pmu.a = measurements.synpmu(:, 6);
meas.pmu.mtrue = measurements.synpmu(:, 7);
meas.pmu.atrue = measurements.synpmu(:, 8);
meas.pmu.Ii = find(meas.pmu.type == 2);
meas.pmu.Iij = find(meas.pmu.type == 1 & meas.pmu.loc > 0);
meas.pmu.Iji = find(meas.pmu.type == 1 & meas.pmu.loc < 0);
meas.pmu.v = find(meas.pmu.type == 3);
% --------------------- Frequency measurements ----------------------------
if isfield(measurements, 'fpmu')
    meas.fpmu.onbus = measurements.fpmu(:, 2);
    meas.fpmu.m = measurements.fpmu(:, 3);
    meas.num.f = numel(meas.fpmu.onbus);
end
% -------------------------------------------------------------------------
[ ~, pmuidx ] = ismember(meas.pmu.onbus, powsys.pmu.onbus); 
meas.pmu.msd = ((powsys.pmu.msd(pmuidx) .* meas.pmu.m)/100);
meas.pmu.asd = (powsys.pmu.asd(pmuidx) * pi / 180);
if isfield(measurements, 'fpmu')
    [ ~, fpmuidx ] = ismember(meas.fpmu.onbus, powsys.pmu.onbus);
	meas.fpmu.fsd = (powsys.pmu.fsd(fpmuidx));
end
% -------------------------------------------------------------------------

% ------------------------ Number of measurements -------------------------
meas.num.scada = size(measurements.scada, 1);
meas.num.pmu = size(measurements.synpmu, 1);

% ----------------------------- SCADA (-s) --------------------------------
meas.num.sPji = numel(meas.scada.pji);
meas.num.sPij = numel(meas.scada.pij);
meas.num.sQji = numel(meas.scada.qji);
meas.num.sQij = numel(meas.scada.qij);
meas.num.sPi = numel(meas.scada.pi);
meas.num.sQi = numel(meas.scada.qi);
meas.num.sIjim = numel(meas.scada.Ijim);
meas.num.sIijm = numel(meas.scada.Iijm);
meas.num.sVm = numel(meas.scada.vm);
meas.num.sSji = numel(meas.scada.sji.p);
meas.num.sSij = numel(meas.scada.sij.p);
meas.num.sSi = numel(meas.scada.si.p);
meas.num.soPji = meas.num.sPji - meas.num.sSji;
meas.num.soQji = meas.num.sQji - meas.num.sSji;
meas.num.soPij = meas.num.sPij - meas.num.sSij;
meas.num.soQij = meas.num.sQij - meas.num.sSij;
meas.num.soPi = meas.num.sPi - meas.num.sSi;
meas.num.soQi = meas.num.sQi - meas.num.sSi;
% -------------------------------------------------------------------------

% ------------------------------ PMU (-p) ---------------------------------
meas.num.pIi = numel(meas.pmu.Ii);
meas.num.pIij = numel(meas.pmu.Iij);
meas.num.pIji = numel(meas.pmu.Iji);
meas.num.pV = numel(meas.pmu.v);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% --------------------------- True values ---------------------------------
meas.Vtrue = measurements.trueVoltage;
% -------------------------------------------------------------------------
end

