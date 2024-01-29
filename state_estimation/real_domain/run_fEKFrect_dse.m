function [ x, x_, converged, info ] = run_fEKFrect_dse(powsys, meas, dsesettings, x_)
info.method = "Extended Kalman Filter with frequency tracking - Rectangular coordinates";
info.paper = [ 'New Kalman Filter Approach Exploiting Frequency ', ...
               'Knowledge for Accurate PMU-based Power System State Estimation' ];
% ------------------------- Variables definition --------------------------
global H Q P_ iAidx jAidx;
reidx = 1:2:2 * powsys.num.bus - 1;
imidx = 2:2:2 * powsys.num.bus;
Trr = 1 / dsesettings.fc;
% -------------------------------------------------------------------------

% ------------ Matrices computed only at the inital run -------------------
if dsesettings.initialStage
    % -------------------- Measurement matrix - H -------------------------
    accI = 0;
    [ rowInj, colInj ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));  
    nnzY = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));
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
    iHIiRe = [ accI + rowInj, accI + rowInj ];
    accI = accI + meas.num.pIi;
    iHIiIm = [ accI + rowInj, accI + rowInj ];
    accI = accI + meas.num.pIi;    
    iHf = (accI + (1:meas.num.f));
    accI = accI + meas.num.f;
    % ---------------------------------------------------------------------
    
    % -------------------------- Column indices ---------------------------
    jHIji = [ 2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji)) - 1; 
               2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji));
               2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)) - 1;
               2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)) ];
    jHIij = [  2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)) - 1; 
               2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.Iij));
               2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.Iij)) - 1;
               2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.Iij)) ];
    jHvRe  = 2 .* meas.pmu.onbus(meas.pmu.v) - 1;
    jHvIm  = 2 .* meas.pmu.onbus(meas.pmu.v);
    jHIiRe = [ (2 * colInj - 1); 2 * colInj ];
    jHIiIm = [ (2 * colInj - 1); 2 * colInj ];
    jHf = (2 * powsys.num.bus + 1) .* ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    
    % ----------------------- Elements' values ----------------------------
    vHIjiRe = [ real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 -imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)));
                 -imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)))];
    vHIjiIm = [ imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
                 imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)));
                 real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)))];
    vHIijRe = [   real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 -imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
                 -imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
    vHIijIm = [ imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
                 imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
    vHv = ones(meas.num.pV, 1);
    vHIiRe = [ real(nnzY), -imag(nnzY) ];
    vHIiIm = [ imag(nnzY),  real(nnzY) ];
    vHf = ones(meas.num.f, 1);
    % ---------------------------------------------------------------------
    if dsesettings.virtual
        [ rowZI, colZI ] = find(powsys.ybus.y(powsys.bus.zi, :));  
        nnzY = nonzeros(powsys.ybus.y(powsys.bus.zi, :));
        % ---------------------- Row indices ------------------------------
        iZI = [ rowZI', rowZI' ];       
        % -----------------------------------------------------------------
        % -------------------- Column indices -----------------------------
        jZI = [ (2 * colZI - 1); 2 * colZI ];
        % -----------------------------------------------------------------
        % -------------------------- Values -------------------------------
        vZI_Re = [ real(nnzY); -imag(nnzY) ];
        vZI_Im = [ imag(nnzY);  real(nnzY) ];
        % -----------------------------------------------------------------
        H = sparse([ iHIjiRe, iHIjiIm, iHIijRe, iHIijIm, iHvRe, iHvIm, iHIiRe', iHIiIm', iHf, accI + iZI, accI + powsys.num.zi + iZI ],...
                   [ jHIji; jHIji; jHIij; jHIij; jHvRe; jHvIm; jHIiRe; jHIiIm; jHf;  jZI; jZI ],...
                   [ vHIjiRe; vHIjiIm; vHIijRe; vHIijIm; vHv; vHv; vHIiRe; vHIiIm; vHf; vZI_Re; vZI_Im ]);
    else
        H = sparse([ iHIjiRe, iHIjiIm, iHIijRe, iHIijIm, iHvRe, iHvIm, iHIiRe, iHIiIm, iHf ],...
                   [ jHIji; jHIji; jHIij; jHIij; jHvRe; jHvIm; jHIiRe; jHIiIm; jHf ],...
                   [ vHIjiRe; vHIjiIm; vHIijRe; vHIijIm; vHv; vHv; vHIiRe; vHIiIm; vHf ]);
    end
    % ---------------------------------------------------------------------

    % -------------------- Process noise covarinace matrix ----------------
    Q = [ 1e-9 .* eye(2 * powsys.num.bus)     zeros(2 * powsys.num.bus, 1) 
            zeros(1, 2 * powsys.num.bus)                 10e-6;
            ];
    % ---------------------------------------------------------------------
    
    % -------------------- Initial prediction covariance ------------------
    P_ = 10 * eye(2 * powsys.num.bus + 1);
    % ---------------------------------------------------------------------
    
    % -------- Row and column indices of the state process Jacobian -------
    iAidx = [ (1:2 * powsys.num.bus), reidx, imidx, reidx, imidx,  2 * powsys.num.bus + 1];
    jAidx = [ (1:2 * powsys.num.bus), imidx, reidx, (2 * powsys.num.bus  + 1)...
               .* ones(1, 2 * powsys.num.bus + 1)];
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

% ------------------------ Measurements vector ----------------------------
z = [
      meas.pmu.m(meas.pmu.Iji) .* cos(meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iji) .* sin(meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* cos(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.Iij) .* sin(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* cos(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.v) .* sin(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Ii) .* cos(meas.pmu.a(meas.pmu.Ii));
      meas.pmu.m(meas.pmu.Ii) .* sin(meas.pmu.a(meas.pmu.Ii));
      meas.fpmu.m - powsys.fn;
      zeros(2 * powsys.num.zi * dsesettings.virtual, 1)
    ];
% -------------------------------------------------------------------------

% -------------- Measurements covariance matrix - R -----------------------
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
vRIi = [ (cos(meas.pmu.a(meas.pmu.Ii)).^2) .* (meas.pmu.msd(meas.pmu.Ii) .^2)...
        + meas.pmu.m(meas.pmu.Ii).^2 .* (sin(meas.pmu.a(meas.pmu.Ii)).^2) .* ...
        (meas.pmu.asd(meas.pmu.Ii) .^2);...
        (sin(meas.pmu.a(meas.pmu.Ii)).^2) .* (meas.pmu.msd(meas.pmu.Ii) .^2)...
        + meas.pmu.m(meas.pmu.Ii).^2 .* (cos(meas.pmu.a(meas.pmu.Ii)).^2) .* ...
        (meas.pmu.asd(meas.pmu.Ii) .^2) 
        ];
vRf = (meas.fpmu.fsd).^2;
vR = [ vRIji; vRIij; vRv; vRIi; vRf ];
% -------------------------------------------------------------------------
if strcmp(dsesettings.mweights(1), "deviceinfo")
    R = [ sparse(1:2 * meas.num.pmu + meas.num.f, 1:2 * meas.num.pmu + meas.num.f, vR),...
        sparse(2 * meas.num.pmu + meas.num.f, 2 * powsys.num.zi * dsesettings.virtual);
        sparse(2 * powsys.num.zi * dsesettings.virtual, 2 * meas.num.pmu + meas.num.f), ...
        sparse(1:2 * powsys.num.zi * dsesettings.virtual,...
        1:2 * powsys.num.zi  * dsesettings.virtual, min(vR)/5 .*...
        ones(2 * powsys.num.zi * dsesettings.virtual, 1))
        ];
else
    R = sparse(1:2 * (meas.num.pmu + powsys.num.zi * dsesettings.virtual + meas.num.f/2),...
               1:2 * (meas.num.pmu + powsys.num.zi * dsesettings.virtual + meas.num.f/2),...
               [ 1e-3 .* ones(2 * meas.num.pmu + meas.num.f, 1); 1e-4 * ones(2 * powsys.num.zi * dsesettings.virtual, 1) ]);
end
% -------------------------------------------------------------------------

% ------------------------ Correction step --------------------------------
K = P_ * H.' * (H * P_ * H.' + R)^-1;
x = x_ + K * (z - H * x_);
P = (eye(2 * powsys.num.bus + 1) - K * H) * P_;
% -------------------------------------------------------------------------

% ------------------------- Prediction step -------------------------------
dA = 2 * pi * x(2 * powsys.num.bus + 1) * Trr;
x_ =  [ reshape([ (x(reidx) * cos(dA) - x(imidx) * sin(dA)).'; (x(reidx) *...
        sin(dA) + x(imidx) * cos(dA)).' ], [], 1); x(2 * powsys.num.bus + 1) ];
% ------------------------ Calculate values of A --------------------------
vA = [   ones(2 * powsys.num.bus, 1); - dA .* ones(powsys.num.bus, 1); 
         dA .* ones(powsys.num.bus, 1);
        [ -x(reidx) .* (2 * pi * Trr)^2 * x(2 * powsys.num.bus + 1) - x(imidx) .* 2 * pi * Trr;
         -x(imidx) .* (2 * pi * Trr)^2 * x(2 * powsys.num.bus + 1) + x(reidx) .* 2 * pi * Trr]; 1 ];
A = sparse(iAidx, jAidx, vA);

P_ = A * P * A.' + Q;
% -------------------------------------------------------------------------
converged = 1;
% -------------------------------------------------------------------------
end

