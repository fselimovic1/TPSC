function [ x, x_, converged, info ] = run_qDKF_dse(powsys, meas, dsesettings, x_)
info.method = "Discrete Kalman Filter with the online assessment of Q";
info.paper = [ 'Performance Assessment of Linear State Estimators ', ...
               'Using Synchrophasor Measurements'...
               ];
% ------------------------- Variables definition --------------------------
converged = 1;
global X H R P_ A;
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
        H = sparse([ iHIjiRe, iHIjiIm, iHIijRe, iHIijIm, iHvRe, iHvIm, iHIiRe', iHIiIm', accI + iZI, accI + powsys.num.zi + iZI ],...
                   [ jHIji; jHIji; jHIij; jHIij; jHvRe; jHvIm; jHIiRe; jHIiIm;  jZI; jZI ],...
                   [ vHIjiRe; vHIjiIm; vHIijRe; vHIijIm; vHv; vHv; vHIiRe; vHIiIm; vZI_Re; vZI_Im ]);
    else
        H = sparse([ iHIjiRe, iHIjiIm, iHIijRe, iHIijIm, iHvRe, iHvIm, iHIiRe, iHIiIm ],...
                   [ jHIji; jHIji; jHIij; jHIij; jHvRe; jHvIm; jHIiRe; jHIiIm ],...
                   [ vHIjiRe; vHIjiIm; vHIijRe; vHIijIm; vHv; vHv; vHIiRe; vHIiIm ]);
    end
    % ---------------------------------------------------------------------

    % ------------------ Measurement noise covarinace matrix --------------
    if strcmp(dsesettings.mweights(1), "deviceinfo")
        % ----------------- Compute values of the elements ----------------
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
        vR = [ vRIji; vRIij; vRv; vRIi ];
        % -----------------------------------------------------------------
        R = [ sparse(1:2 * meas.num.pmu, 1:2 * meas.num.pmu, vR), ...
              sparse(2 * meas.num.pmu, 2 * powsys.num.zi * dsesettings.virtual);
              sparse(2 * powsys.num.zi * dsesettings.virtual, 2 * meas.num.pmu), ...
              sparse(1:2 * powsys.num.zi * dsesettings.virtual,...
              1:2 * powsys.num.zi  * dsesettings.virtual, min(vR)/5 .*...
              ones(2 * powsys.num.zi * dsesettings.virtual, 1))
            ];
        % -----------------------------------------------------------------
    else
        R = sparse(1:2 * (meas.num.pmu + powsys.num.zi * dsesettings.virtual),...
                   1:2 * (meas.num.pmu + powsys.num.zi * dsesettings.virtual),...
                   [ ones(2 * meas.num.pmu, 1); 5 * ones(2 * powsys.num.zi * dsesettings.virtual, 1) ]);
    end
    % -------------------- State Transition Matrix ------------------------
    A = sparse(1:2 * powsys.num.bus, 1:2 * powsys.num.bus, ones(2 * powsys.num.bus, 1));
    % ---------------------------------------------------------------------
    
    % -------------------- Initial prediction covariance ------------------
    P_ = 10 * eye(2 * powsys.num.bus);
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
      zeros(2 * powsys.num.zi * dsesettings.virtual, 1)
    ];
% -------------------------------------------------------------------------


% -------------------- Update measurements covariance matrix --------------
if dsesettings.timevariantR && strcmp(dsesettings.mweights(1), "deviceinfo")
    % ----------------- Compute values of the elements --------------------
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
    vR = [ vRIji; vRIij; vRv; vRIi ];
    % ---------------------------------------------------------------------
    R = [ sparse(1:2 * meas.num.pmu, 1:2 * meas.num.pmu, vR), ...
          sparse(2 * meas.num.pmu, 2 * powsys.num.zi * dsesettings.virtual);
          sparse(2 * powsys.num.zi * dsesettings.virtual, 2 * meas.num.pmu), ...
          sparse(1:2 * powsys.num.zi * dsesettings.virtual,...
          1:2 * powsys.num.zi  * dsesettings.virtual, min(vR)/5 .*...
          ones(2 * powsys.num.zi * dsesettings.virtual, 1))
            ];
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

% ------------------- Calculate process covariance matrix -----------------
if dsesettings.tstep == 1
    Q = A .* 10^-9;
elseif dsesettings.tstep <= dsesettings.NQ
    Q = diag(var(X(:, 1:dsesettings.tstep - 1), 0, 2));
else
    Q = diag(var(X(:, dsesettings.tstep - dsesettings.NQ:dsesettings.tstep), 0, 2));
end
% -------------------------------------------------------------------------

% ----------------------- KF COMPUTING PART -------------------------------
% -------------------------------------------------------------------------

% ------------------------ Correction step --------------------------------
K = P_ * H.' * (H * P_ * H.' + R)^-1;
x = x_ + K * (z - H * x_);
P = (A - K * H) * P_;
% -------------------------------------------------------------------------

% ------------------------- Prediction step -------------------------------
x_ =  x;
P_ = P + Q;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
end

