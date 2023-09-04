function [ X, X_, converged, info ] = run_fEKFrect_rtse(powsys, meas, rtsesettings, X_)
info.method = "Extended Kalman Filter with frequency tracking - rectangular coordinates";
info.paper = [ 'New Kalman Filter Approach Exploiting Frequency ', ...
     'Knowledge for Accurate PMU-based Power System State Estimation' ];
% ------------------------- Variables definition --------------------------
global H Q P_ iAidx jAidx iRidx jRidx;
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
    jHibrO = [ 2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.IijO)) - 1; 
               2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.IijO));
               2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.IijO)) - 1;
               2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.IijO)) ];
    jHibr = [  2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)) - 1; 
               2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.Iij));
               2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.Iij)) - 1;
               2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.Iij)) ];
    jHvRe  = 2 .* meas.pmu.onbus(meas.pmu.v) - 1;
    jHvIm  = 2 .* meas.pmu.onbus(meas.pmu.v);
    jHf = (2 * powsys.num.bus + 1) .* ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    
    % ----------------------- Elements' values ----------------------------
    vHibrOre = [ real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.IijO)));
                 -imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.IijO)));
                 real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.IijO)));
                 -imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.IijO)))];
    vHibrOim = [ imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.IijO)));
                 real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.IijO)));
                 imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.IijO)));
                 real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.IijO)))];
    vHibrRe = [   real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 -imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
                 -imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
    vHibrIm = [ imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
    vHv = ones(meas.num.pV, 1);
    vHf = ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    H = sparse([ iHibrOre, iHibrOim, iHibrRe, iHibrIm, iHvRe, iHvIm, iHf ],...
               [ jHibrO; jHibrO; jHibr; jHibr; jHvRe; jHvIm; jHf ],...
               [ vHibrOre; vHibrOim; vHibrRe; vHibrIm; vHv; vHv; vHf ]);
    % ---------------------------------------------------------------------
    
%     % ------------------ Measurement noise covarinace matrix --------------
    % --------------------------- Row indices -----------------------------
    accI = meas.num.pIijO;
    iRIbrO = [ (1:meas.num.pIijO), (1:meas.num.pIijO), ...
               accI + (1:meas.num.pIijO), accI + (1:meas.num.pIijO)];      
    accI = accI + meas.num.pIijO;
    iRIbr = [ accI + (1:meas.num.pIij), accI + (1:meas.num.pIij),...
              accI + meas.num.pIij + (1:meas.num.pIij), accI + meas.num.pIij + (1:meas.num.pIij) ];
    accI = accI + 2 * meas.num.pIij;
    iRv = [ accI + (1:meas.num.pV), accI + (1:meas.num.pV),...
            accI + meas.num.pV + (1:meas.num.pV), accI + meas.num.pV + (1:meas.num.pV) ];
    accI = accI + 2 * meas.num.pV;
    iRf =  accI + (1:meas.num.f);
    iRidx = [ iRIbrO, iRIbr, iRv, iRf ];
    % ---------------------------------------------------------------------
    
    % ------------------------- Column indices ----------------------------
    accJ = meas.num.pIijO;
    jRIbrO = [ 1:meas.num.pIijO, accJ + (1:meas.num.pIijO), ...
               accJ + (1:meas.num.pIijO), 1:meas.num.pIijO ];
    accJ = accJ + meas.num.pIijO;
    jRIbr = [ (accJ + (1:meas.num.pIij)), accJ + meas.num.pIij + (1:meas.num.pIij), ...
               accJ + meas.num.pIij + (1:meas.num.pIij), (accJ + (1:meas.num.pIij)) ];
    accJ = accJ + 2 * meas.num.pIij;   
    jRv = [ (accJ + (1:meas.num.pV)), accJ + meas.num.pV + (1:meas.num.pV), ...
             accJ + meas.num.pV + (1:meas.num.pV), (accJ + (1:meas.num.pV)) ];
    accJ = accJ + 2 * meas.num.pV;
    jRf = accJ + (1:meas.num.f);
    jRidx = [ jRIbrO, jRIbr, jRv, jRf ];
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    % -------------------- Process noise covarinace matrix ----------------
    Q = 10^-9.* eye(2 * powsys.num.bus + 1);
    % ---------------------------------------------------------------------
    
    % -------------------- Initial prediction covariance ------------------
    P_ = 10 * eye(2 * powsys.num.bus + 1);
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
      meas.pmu.m(meas.pmu.IijO) .* cos(meas.pmu.a(meas.pmu.IijO));
      meas.pmu.m(meas.pmu.IijO) .* sin(meas.pmu.a(meas.pmu.IijO));
      meas.pmu.m(meas.pmu.Iij) .* cos(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.Iij) .* sin(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* cos(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.v) .* sin(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* cos(meas.pmu.a(meas.pmu.Iinj));
      meas.pmu.m(meas.pmu.Iinj) .* sin(meas.pmu.a(meas.pmu.Iinj));
      meas.fpmu.m - powsys.fn;
    ];
% -------------------------------------------------------------------------

% -------------- Measurements covariance matrix - R -----------------------
% ----------------- Compute values of the elements ------------------------
vRIbrO = [ (cos(meas.pmu.a(meas.pmu.IijO)).^2) .* (meas.pmu.msd(meas.pmu.IijO) .^2)...
            + meas.pmu.m(meas.pmu.IijO).^2 .* (sin(meas.pmu.a(meas.pmu.IijO)).^2) .* ...
            (meas.pmu.asd(meas.pmu.IijO) .^2);...
            cos(meas.pmu.a(meas.pmu.IijO)) .* sin(meas.pmu.a(meas.pmu.IijO)) .* ...
            (meas.pmu.msd(meas.pmu.IijO) .^2) - meas.pmu.m(meas.pmu.IijO).^2 .* ...
            sin(meas.pmu.a(meas.pmu.IijO)) .* cos(meas.pmu.a(meas.pmu.IijO)) ...
            .* (meas.pmu.asd(meas.pmu.IijO) .^2);
            (sin(meas.pmu.a(meas.pmu.IijO)).^2) .* (meas.pmu.msd(meas.pmu.IijO) .^2)...
            + meas.pmu.m(meas.pmu.IijO).^2 .* (cos(meas.pmu.a(meas.pmu.IijO)).^2) .* ...
            (meas.pmu.asd(meas.pmu.IijO) .^2);
            cos(meas.pmu.a(meas.pmu.IijO)) .* sin(meas.pmu.a(meas.pmu.IijO)) .* ...
            (meas.pmu.msd(meas.pmu.IijO) .^2) - meas.pmu.m(meas.pmu.IijO).^2 .* ...
            sin(meas.pmu.a(meas.pmu.IijO)) .* cos(meas.pmu.a(meas.pmu.IijO)) ...
            .* (meas.pmu.asd(meas.pmu.IijO) .^2);
             ];      
vRIbr = [ (cos(meas.pmu.a(meas.pmu.Iij)).^2) .* (meas.pmu.msd(meas.pmu.Iij) .^2)...
            + meas.pmu.m(meas.pmu.Iij).^2 .* (sin(meas.pmu.a(meas.pmu.Iij)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);...
            cos(meas.pmu.a(meas.pmu.Iij)) .* sin(meas.pmu.a(meas.pmu.Iij)) .* ...
            (meas.pmu.msd(meas.pmu.Iij) .^2) - meas.pmu.m(meas.pmu.Iij).^2 .* ...
            sin(meas.pmu.a(meas.pmu.Iij)) .* cos(meas.pmu.a(meas.pmu.Iij)) ...
            .* (meas.pmu.asd(meas.pmu.Iij) .^2);
            (sin(meas.pmu.a(meas.pmu.Iij)).^2) .* (meas.pmu.msd(meas.pmu.Iij) .^2)...
            + meas.pmu.m(meas.pmu.Iij).^2 .* (cos(meas.pmu.a(meas.pmu.Iij)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);
            cos(meas.pmu.a(meas.pmu.Iij)) .* sin(meas.pmu.a(meas.pmu.Iij)) .* ...
            (meas.pmu.msd(meas.pmu.Iij) .^2) - meas.pmu.m(meas.pmu.Iij).^2 .* ...
            sin(meas.pmu.a(meas.pmu.Iij)) .* cos(meas.pmu.a(meas.pmu.Iij)) ...
            .* (meas.pmu.asd(meas.pmu.Iij) .^2);
             ]; 
vRv = [ (cos(meas.pmu.a(meas.pmu.v)).^2) .* (meas.pmu.msd(meas.pmu.v) .^2)...
            + meas.pmu.m(meas.pmu.v).^2 .* (sin(meas.pmu.a(meas.pmu.v)).^2) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);...
            cos(meas.pmu.a(meas.pmu.v)) .* sin(meas.pmu.a(meas.pmu.v)) .* ...
            (meas.pmu.msd(meas.pmu.v) .^2) - meas.pmu.m(meas.pmu.v).^2 .* ...
            sin(meas.pmu.a(meas.pmu.v)) .* cos(meas.pmu.a(meas.pmu.v)) ...
            .* (meas.pmu.asd(meas.pmu.v) .^2);
            (sin(meas.pmu.a(meas.pmu.v)).^2) .* (meas.pmu.msd(meas.pmu.v) .^2)...
            + meas.pmu.m(meas.pmu.v).^2 .* (cos(meas.pmu.a(meas.pmu.v)).^2) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);
            cos(meas.pmu.a(meas.pmu.v)) .* sin(meas.pmu.a(meas.pmu.v)) .* ...
            (meas.pmu.msd(meas.pmu.v) .^2) - meas.pmu.m(meas.pmu.v).^2 .* ...
            sin(meas.pmu.a(meas.pmu.v)) .* cos(meas.pmu.a(meas.pmu.v)) ...
            .* (meas.pmu.asd(meas.pmu.v) .^2);
             ];    
vRf = (meas.fpmu.fsd).^2;
vR = [ vRIbrO; vRIbr; vRv; vRf ];
% -------------------------------------------------------------------------
R = sparse(iRidx, jRidx, vR);
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

