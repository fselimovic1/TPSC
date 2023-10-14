function [ X, X_, converged, info ] = run_fEKFrect_dse(powsys, meas, dsesettings, X_)
info.method = "Extended Kalman Filter with fRequency tracking - Rectangular coordinates";
info.paper = [ 'New Kalman Filter Approach Exploiting FRequency ', ...
     'Knowledge for Accurate PMU-based Power System State EstImation' ];
% ------------------------- Variables definition --------------------------
global H Q P_ iAidx jAidx iRidx jRidx;
Reidx = 1:2:2 * powsys.num.bus - 1;
Imidx = 2:2:2 * powsys.num.bus;
Trr = 1 / dsesettings.fc;
% -------------------------------------------------------------------------

% ------------ Matrices computed only at the inital run -------------------
if dsesettings.initialStage
    % -------------------- Measurement matrix - H -------------------------
    accI = 0;
    % ------------------------ Row indices --------------------------------
    iHIjiRe = [ 1:meas.num.pIji, 1:meas.num.pIji, 1:meas.num.pIji, 1:meas.num.pIji ];
    accI = accI + meas.num.pIji;
    iHIjiIm = [ (accI + (1:meas.num.pIji)), (accI + (1:meas.num.pIji)),...
                (accI + (1:meas.num.pIji)), (accI + (1:meas.num.pIji)) ];
    accI = accI + meas.num.pIji; 
    iHIijRe = [ (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)),...
                (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)) ];
    accI = accI + meas.num.pIij;         
    iHIijIm = [ (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)),...
                (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)) ]; 
    accI = accI + meas.num.pIij; 
    iHvRe = (accI + (1:meas.num.pV));
    accI = accI + meas.num.pV; 
    iHvIm = (accI + (1:meas.num.pV));
    accI = accI + meas.num.pV;
    iHf = (accI + (1:meas.num.f));
    % ---------------------------------------------------------------------
    
    % -------------------------- Column indices ---------------------------
    jHIji = [ 2 .* powsys.ijanch.i(-meas.pmu.loc(meas.pmu.Iji)) - 1; 
               2 .* powsys.ijanch.i(-meas.pmu.loc(meas.pmu.Iji));
               2 .* powsys.ijanch.j(-meas.pmu.loc(meas.pmu.Iji)) - 1;
               2 .* powsys.ijanch.j(-meas.pmu.loc(meas.pmu.Iji)) ];
    jHIij = [  2 .* powsys.ijanch.i(meas.pmu.loc(meas.pmu.Iij)) - 1; 
               2 .* powsys.ijanch.i(meas.pmu.loc(meas.pmu.Iij));
               2 .* powsys.ijanch.j(meas.pmu.loc(meas.pmu.Iij)) - 1;
               2 .* powsys.ijanch.j(meas.pmu.loc(meas.pmu.Iij)) ];
    jHvRe  = 2 .* meas.pmu.onbus(meas.pmu.v) - 1;
    jHvIm  = 2 .* meas.pmu.onbus(meas.pmu.v);
    jHf = (2 * powsys.num.bus + 1) .* ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    
    % ----------------------- Elements' values ----------------------------
    vHIjiRe = [ Real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 -Imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 Real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)));
                 -Imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)))];
    vHIjiIm = [ Imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 Real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 Imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)));
                 Real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)))];
    vHIijRe = [   Real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 -Imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 Real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
                 -Imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
    vHIijIm = [ Imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 Real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 Imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
                 Real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
    vHv = ones(meas.num.pV, 1);
    vHf = ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    H = sparse([ iHIjiRe, iHIjiIm, iHIijRe, iHIijIm, iHvRe, iHvIm, iHf ],...
               [ jHIji; jHIji; jHIij; jHIij; jHvRe; jHvIm; jHf ],...
               [ vHIjiRe; vHIjiIm; vHIijRe; vHIijIm; vHv; vHv; vHf ]);
    % ---------------------------------------------------------------------
    
    % ------------------ MeasuRement noise covarinace matrix --------------
    % ----------------------------- Row indices ---------------------------
    accI = meas.num.pIji;
    iRIji = [ (1:meas.num.pIji), accI + (1:meas.num.pIji)]; 
    accI = accI + meas.num.pIji;
    iRIij = [ accI + (1:meas.num.pIij), accI + meas.num.pIij + (1:meas.num.pIij)];
    accI = accI + 2 * meas.num.pIij;
    iRv = [ accI + (1:meas.num.pV), accI + meas.num.pV + (1:meas.num.pV)];
    accI = accI + 2 * meas.num.pV;
    iRf =  accI + (1:meas.num.f);
    iRidx = [ iRIji, iRIij, iRv, iRf ];
    % ---------------------------------------------------------------------

    % ---------------------------- Column indices -------------------------
    accJ = meas.num.pIji;
    jRIji = [ 1:meas.num.pIji, accJ + (1:meas.num.pIji) ];
    accJ = accJ + meas.num.pIji;
    jRIij = [ (accJ + (1:meas.num.pIij)), accJ + meas.num.pIij + (1:meas.num.pIij)];
    accJ = accJ + 2 * meas.num.pIij;   
    jRv = [ (accJ + (1:meas.num.pV)), accJ + meas.num.pV + (1:meas.num.pV)];
    accJ = accJ + 2 * meas.num.pV;
    jRf = accJ + (1:meas.num.f);
    jRidx = [ jRIji, jRIij, jRv, jRf ];
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------

    % -------------------- Process noise covarinace matrix ----------------
    Q = 10^-4.* eye(2 * powsys.num.bus + 1);
    % ---------------------------------------------------------------------
    
    % -------------------- Initial pRediction covariance ------------------
    P_ = 10 * eye(2 * powsys.num.bus + 1);
    % ---------------------------------------------------------------------
    
    % -------- Row and column indices of the state process Jacobian -------
    iAidx = [ (1:2 * powsys.num.bus), (1:2 * powsys.num.bus), Reidx, Imidx,  2 * powsys.num.bus + 1];
    jAidx = [ (1:2 * powsys.num.bus), Reshape([Imidx; ...
               Reidx], 1, []), (2 * powsys.num.bus  + 1)...
               .* ones(1, 2 * powsys.num.bus + 1)];
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

% ------------------------ MeasuRements vector ----------------------------
z = [
      meas.pmu.m(meas.pmu.Iji) .* cos(meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iji) .* sin(meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* cos(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.Iij) .* sin(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* cos(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.v) .* sin(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* cos(meas.pmu.a(meas.pmu.Iinj));
      meas.pmu.m(meas.pmu.Iinj) .* sin(meas.pmu.a(meas.pmu.Iinj));
      meas.fpmu.m - powsys.fn;
    ];
% -------------------------------------------------------------------------

% -------------- MeasuRements covariance matrix - R -----------------------
% ----------------- Compute values of the elements ------------------------
vRIji = [ (cos(meas.pmu.a(meas.pmu.Iji)).^2) .* (meas.pmu.msd(meas.pmu.Iji) .^2)...
            + meas.pmu.m(meas.pmu.Iji).^2 .* (sin(meas.pmu.a(meas.pmu.Iji)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);...
            (sin(meas.pmu.a(meas.pmu.Iji)).^2) .* (meas.pmu.msd(meas.pmu.Iji) .^2)...
            + meas.pmu.m(meas.pmu.Iji).^2 .* (cos(meas.pmu.a(meas.pmu.Iji)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);
             ];      
vRIij = [ (cos(meas.pmu.a(meas.pmu.Iij)).^2) .* (meas.pmu.msd(meas.pmu.Iij) .^2)...
            + meas.pmu.m(meas.pmu.Iij).^2 .* (sin(meas.pmu.a(meas.pmu.Iij)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);...
            (sin(meas.pmu.a(meas.pmu.Iij)).^2) .* (meas.pmu.msd(meas.pmu.Iij) .^2)...
            + meas.pmu.m(meas.pmu.Iij).^2 .* (cos(meas.pmu.a(meas.pmu.Iij)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);
             ]; 
vRv = [ (cos(meas.pmu.a(meas.pmu.v)).^2) .* (meas.pmu.msd(meas.pmu.v) .^2)...
            + meas.pmu.m(meas.pmu.v).^2 .* (sin(meas.pmu.a(meas.pmu.v)).^2) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);...
            (sin(meas.pmu.a(meas.pmu.v)).^2) .* (meas.pmu.msd(meas.pmu.v) .^2)...
            + meas.pmu.m(meas.pmu.v).^2 .* (cos(meas.pmu.a(meas.pmu.v)).^2) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);
             ];    
vRf = (meas.fpmu.fsd).^2;
vR = [ vRIji; vRIij; vRv; vRf ];
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
X_ =  [ Reshape([ (X(Reidx) * cos(dA) - X(Imidx) * sin(dA)).'; (X(Reidx) *...
        sin(dA) + X(Imidx) * cos(dA)).' ], [], 1); X(2 * powsys.num.bus + 1) ];
% ------------------------ Calculate values of A --------------------------
vA = [   ones(2 * powsys.num.bus, 1); - dA .* ones(2 * powsys.num.bus, 1); ...
       [ -X(Reidx) .* (2 * pi * Trr)^2 * X(2 * powsys.num.bus + 1) - X(Imidx) .* 2 * pi * Trr;
         -X(Imidx) .* (2 * pi * Trr)^2 * X(2 * powsys.num.bus + 1) + X(Reidx) .* 2 * pi * Trr]; 1 ];
A = sparse(iAidx, jAidx, vA);

P_ = A * P * A.' + Q;
% -------------------------------------------------------------------------
converged = 1;
% -------------------------------------------------------------------------
end

