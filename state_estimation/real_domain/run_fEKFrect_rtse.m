function [ X, X_, converged, info ] = run_fEKFrect_rtse(powsys, meas, rtsesettings, X_)
info.method = "Extended Kalman Filter with frequency tracking - rectangular coordinates";
info.paper = [ 'New Kalman Filter Approach Exploiting Frequency ', ...
     'Knowledge for Accurate PMU-based Power System State Estimation' ];
% ------------------------- Variables definition --------------------------
global H R Q P_ iAidx jAidx;
reidx = 1:2:2 * powsys.num.bus - 1;
imidx = 2:2:2 * powsys.num.bus;
Trr = 1 / rtsesettings.fc;
% -------------------------------------------------------------------------

% ------------ Matrices computed only at the inital run -------------------
if rtsesettings.initialStage
    % -------------------- Measurement matrix - H -------------------------
    accI = 0;
    % ------------------------ Row indices --------------------------------
    iHibrOre = [ 1:meas.num.pIijO, 1:meas.num.pIijO, 1:meas.num.pIijO, 1:meas.num.pIijO ];
    accI = accI + meas.num.pIijO;
    iHibrOim = [ (accI + (1:meas.num.pIijO)), (accI + (1:meas.num.pIijO)),...
                (accI + (1:meas.num.pIijO)), (accI + (1:meas.num.pIijO)) ];
    accI = accI + meas.num.pIijO; 
    iHibrRe = [ (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)),...
                (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)) ];
    accI = accI + meas.num.pIij;         
    iHibrIm = [ (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)),...
                (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)) ]; 
    accI = accI + meas.num.pIij; 
    iHvRe = (accI + (1:meas.num.pV));
    accI = accI + meas.num.pV; 
    iHvIm = (accI + (1:meas.num.pV));
    accI = accI + meas.num.pV;
    iHf = (accI + (1:meas.num.f));
    % ---------------------------------------------------------------------
    
    % -------------------------- Column indices ---------------------------
    jHibrO = [ 2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.ibranchO)) - 1; 
               2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.ibranchO));
               2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.ibranchO)) - 1;
               2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.ibranchO)) ];
    jHibr = [  2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.ibranch)) - 1; 
               2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.ibranch));
               2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.ibranch)) - 1;
               2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.ibranch)) ];
    jHvRe  = 2 .* meas.pmu.onbus(meas.pmu.vnode) - 1;
    jHvIm  = 2 .* meas.pmu.onbus(meas.pmu.vnode);
    jHf = (2 * powsys.num.bus + 1) .* ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    
    % ----------------------- Elements' values ----------------------------
    vHibrOre = [ real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchO)));
                 -imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchO)));
                 real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchO)));
                 -imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchO)))];
    vHibrOim = [ imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchO)));
                 real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchO)));
                 imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchO)));
                 real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchO)))];
    vHibrRe = [   real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 -imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)));
                 -imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)))];
    vHibrIm = [ imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)))];
    vHv = ones(meas.num.pV, 1);
    vHf = ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    H = sparse([ iHibrOre, iHibrOim, iHibrRe, iHibrIm, iHvRe, iHvIm, iHf ],...
               [ jHibrO; jHibrO; jHibr; jHibr; jHvRe; jHvIm; jHf ],...
               [ vHibrOre; vHibrOim; vHibrRe; vHibrIm; vHv; vHv; vHf ]);
    % ---------------------------------------------------------------------
    
    % ------------------ Measurement noise covarinace matrix --------------
    R = 10^-9 .* eye(2 * meas.num.pmu + meas.num.f);
    % ---------------------------------------------------------------------
    
    % -------------------- Process noise covarinace matrix ----------------
    Q = 10^-9 .* eye(2 * powsys.num.bus + 1);
    % ---------------------------------------------------------------------
    
    % -------------------- Initial prediction covariance ------------------
    P_ = eye(2 * powsys.num.bus + 1);
    % ---------------------------------------------------------------------
    
    % -------- Row and column indices of the state process Jacobian -------
    iAidx = [ (1:2 * powsys.num.bus), (1:2 * powsys.num.bus), reidx, imidx,  2 * powsys.num.bus + 1];
    jAidx = [ (1:2 * powsys.num.bus), reshape([imidx; ...
               reidx], 1, []), (2 * powsys.num.bus  + 1)...
               .* ones(1, 2 * powsys.num.bus + 1)];
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

% ------------------------ Measurements vector ----------------------------
z = [
      meas.pmu.m(meas.pmu.ibranchO) .* cos(meas.pmu.a(meas.pmu.ibranchO));
      meas.pmu.m(meas.pmu.ibranchO) .* sin(meas.pmu.a(meas.pmu.ibranchO));
      meas.pmu.m(meas.pmu.ibranch) .* cos(meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.ibranch) .* sin(meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* cos(meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.vnode) .* sin(meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.inj) .* cos(meas.pmu.a(meas.pmu.inj));
      meas.pmu.m(meas.pmu.inj) .* sin(meas.pmu.a(meas.pmu.inj));
      meas.fpmu.m - powsys.fn;
    ];
% -------------------------------------------------------------------------

% ------------------------ Correction step --------------------------------
K = P_ * H.' * (H * P_ * H.' + R)^-1;
X = X_ + K * (z - H * X_);
P = (eye(2 * powsys.num.bus + 1) - K * H) * P_;
% -------------------------------------------------------------------------

% ------------------------- Prediction step -------------------------------
dA = 2 * pi * X(2 * powsys.num.bus + 1) * Trr;
X_ =  [ reshape([ (X(reidx) * cos(dA) - X(imidx) * sin(dA)).'; (X(reidx) *...
        sin(dA) + X(imidx) * cos(dA)).' ], [], 1); X(2 * powsys.num.bus + 1) ];
% ------------------------ Calculate values of A --------------------------
vA = [ ones(2 * powsys.num.bus, 1); - dA .* ones(2 * powsys.num.bus, 1); ...
       [ -X(reidx) .* (2 * pi * Trr)^2 * X(2 * powsys.num.bus + 1) - X(imidx) .* 2 * pi * Trr;
         -X(imidx) .* (2 * pi * Trr)^2 * X(2 * powsys.num.bus + 1) + X(reidx) .* 2 * pi * Trr]; 1 ];
A = sparse(iAidx, jAidx, vA);

P_ = A * P * A.' + Q;
% -------------------------------------------------------------------------
converged = 1;
% -------------------------------------------------------------------------
end

