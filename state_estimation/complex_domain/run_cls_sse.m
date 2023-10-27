function [ Vc, iter, converged, info ] = run_cls_sse(powsys, meas, ssesettings)
info.method = 'Linear State Estimation (PMU only) in Complex Variables';
% --------------------- Setting additional variables ----------------------
iter = 1;
converged = 1;
% -------------------------------------------------------------------------

% ----------------------- Vector of measurement values --------------------
z = [ meas.pmu.m(meas.pmu.Iji) .* exp(1i .* meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Ii) .* exp(1i .* meas.pmu.a(meas.pmu.Ii));
      zeros(powsys.num.zi * ssesettings.virtual, 1)
    ];
% -------------------------------------------------------------------------

% ------------------------- Construct H matrix --------------------------
[ rowi, coli ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));  

% ---------------------------- Row indices --------------------------------
iHIji= [ 1:meas.num.pIji, 1:meas.num.pIji];
iHIij = [ (meas.num.pIji + (1:meas.num.pIij)), (meas.num.pIji + (1:meas.num.pIij)) ];
iHV = (meas.num.pIji+ meas.num.pIij + (1:meas.num.pV));
iHIi = meas.num.pIji+ meas.num.pIij + meas.num.pV + rowi;
accI = max([iHIji, iHIij, iHV, iHIi]);
% -------------------------------------------------------------------------

% ------------------------ Column indices ---------------------------------
jHIji = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)); ];
jHIij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.Iij));
                 ];
jHV = meas.pmu.onbus(meas.pmu.v);
jHIi = coli;
% -------------------------------------------------------------------------

% ------------------------- Values of elements ----------------------------
vHIji= [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)); ];
vHIij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); 
    ];
vHV = ones(meas.num.pV, 1);
vHIi = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));
% -------------------------------------------------------------------------

% ------------------- Virtual measurements (optional) -----------------
if ssesettings.virtual
    [ iZI, jZI, vZI ] = find(powsys.ybus.y(powsys.bus.zi, :));
    H = sparse(...
                 [ iHIji, iHIij, iHV, iHIi, accI + iZI'  ],... 
                 [ jHIji; jHIij; jHV; jHIi; jZI  ], ...
                 [ vHIji; vHIij; vHV; vHIi; vZI  ], ...
                 meas.num.pmu + powsys.num.zi, powsys.num.bus ...
             );
    % --------------------- Construct R matrix -------------------------
    W = sparse(1:meas.num.pmu + powsys.num.zi, 1:meas.num.pmu + powsys.num.zi, ...
        [1e-6 .* ones(meas.num.pmu, 1);  5 * (1e-6) .* ones(powsys.num.zi, 1)]);
    % -----------------------------------------------------------------          
else
    H = sparse(...
                 [ iHIji, iHIij, iHV, iHIi  ],... 
                 [ jHIji; jHIij; jHV; jHIi  ], ...
                 [ vHIji; vHIij; vHV; vHIi  ], ...
                 meas.num.pmu, powsys.num.bus ...
             );
    % --------------------- Construct R matrix ------------------------
    W = sparse(1:meas.num.pmu, 1:meas.num.pmu, 1e-6 .* ones(meas.num.pmu, 1));
    % -----------------------------------------------------------------
end
% -------------------------------------------------------------------------


% --------------------------- Solve the LS problem ------------------------
    x = (H' * W * H) \ (H' * W * z);
% -------------------------------------------------------------------------

% ----------------------- Estimator results -------------------------------
Vc = x;
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end