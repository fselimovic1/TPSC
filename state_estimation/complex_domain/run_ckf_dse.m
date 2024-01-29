function [ x, x_, converged, info ] = run_ckf_dse(powsys, meas, dsesettings, x_)
info.method = "Kalman Filter in Complex Variables";
info.paper = [ ];

% ------------------------- Variables definition --------------------------
global X H R P_ Im
converged = 1;
% -------------------------------------------------------------------------

% -------------------------- PREPARATION PHASE ----------------------------
% -------------------------------------------------------------------------
if dsesettings.initialStage
    % ------------------------- Construct H matrix ------------------------
    [ rowi, coli ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));  

    % ---------------------------- Row indices ----------------------------
    accI = 0;
    iHIji= [ 1:meas.num.pIji, 1:meas.num.pIji];
    accI = accI + meas.num.pIji;
    iHIij = [ accI + (1:meas.num.pIij), accI + (1:meas.num.pIij) ];
    accI = accI + meas.num.pIij;
    iHV = accI + (1:meas.num.pV);
    accI = accI + meas.num.pV;
    iHIi = accI + rowi;
    if ~isempty(rowi)
        accI = accI + max(rowi);
    end
    % ---------------------------------------------------------------------

    % ------------------------ Column indices -----------------------------
    jHIji= [ powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji)); 
             powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)); ];
    jHIij =  [ powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
               powsys.branch.j(meas.pmu.loc(meas.pmu.Iij));
                     ];
    jHV = meas.pmu.onbus(meas.pmu.v);
    jHIi = coli;
    % ---------------------------------------------------------------------

    % ------------------------- Values of elements ------------------------
    vHIji= [ powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji));...
             powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)); ];
    vHIij =  [ powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
               powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); 
        ];
    vHV = ones(meas.num.pV, 1);
    vHIi = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));
    % ---------------------------------------------------------------------
    
    % ------------------- Virtual measurements (optional) -----------------
    if dsesettings.virtual
        [ iZI, jZI, vZI ] = find(powsys.ybus.y(powsys.bus.zi, :));
        H = sparse(...
                     [ iHIji, iHIij, iHV, iHIi, accI + iZI'  ],... 
                     [ jHIji; jHIij; jHV; jHIi; jZI  ], ...
                     [ vHIji; vHIij; vHV; vHIi; vZI  ], ...
                     meas.num.pmu + powsys.num.zi, powsys.num.bus ...
                 );
        % --------------------- Construct R matrix -------------------------
        R = sparse(1:meas.num.pmu + powsys.num.zi, 1:meas.num.pmu + powsys.num.zi, ...
            [1e-4 .* ones(meas.num.pmu, 1);  (1e-5) .* ones(powsys.num.zi, 1)]);
        % -----------------------------------------------------------------          
    else
        H = sparse(...
                     [ iHIji, iHIij, iHV, iHIi  ],... 
                     [ jHIji; jHIij; jHV; jHIi  ], ...
                     [ vHIji; vHIij; vHV; vHIi  ], ...
                     meas.num.pmu, powsys.num.bus ...
                 );
        % --------------------- Construct R matrix ------------------------
        R = sparse(1:meas.num.pmu, 1:meas.num.pmu, 1e-6 .* ones(meas.num.pmu, 1));
        % -----------------------------------------------------------------
    end
    % ---------------------------------------------------------------------
    % -------------------- Initial prediction covariance ------------------
    P_ = 10 * eye(powsys.num.bus);
    % ---------------------------------------------------------------------
    Im = sparse(1:powsys.num.bus, 1:powsys.num.bus, ones(powsys.num.bus, 1));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% --------------------------- REAL TIME PHASE -----------------------------
% -------------------------------------------------------------------------

% ----------------------- Vector of measurement values --------------------
z = [ meas.pmu.m(meas.pmu.Iji) .* exp(1i .* meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Ii) .* exp(1i .* meas.pmu.a(meas.pmu.Ii));
      zeros(powsys.num.zi * dsesettings.virtual, 1);
    ];
% -------------------------------------------------------------------------

% ------------------- Calculate process covariance matrix -----------------
if dsesettings.tstep < 3
    Q = Im .* 10^-9;
elseif dsesettings.tstep <= dsesettings.NQ
    y = X(:, dsesettings.tstep - 1:-1:2) - X(:, 1);
    Q = diag(var(y, 0, 2));
else
    y = X(:, dsesettings.tstep - 1:-1:dsesettings.tstep - dsesettings.NQ + 1)...
        - X(:, dsesettings.tstep - dsesettings.NQ);
    Q = diag(var(y, 0, 2));
end
% -------------------------------------------------------------------------

% ------------------------ Correction step --------------------------------
K = P_ * H.' * (H * P_ * H.' + R)^-1;
x = x_ + K * (z - H * x_);
P = (Im - K * H) * P_;
% -------------------------------------------------------------------------

% ------------------------- Prediction step -------------------------------
x_ =  x;
P_ = P + Q;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
end

