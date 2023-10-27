function [ Vc, iter, converged, info ] = run_lcec_sse(powsys, meas)
info.method = 'Linear Equality Constrained State Estimation (PMU only) in Complex Variables';
% --------------------- Setting additional variables ----------------------
iter = 1;
converged = 1;
% -------------------------------------------------------------------------

% ----------------------- Vector of measurement values --------------------
z = [ meas.pmu.m(meas.pmu.Iji) .* exp(1i .* meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Ii) .* exp(1i .* meas.pmu.a(meas.pmu.Ii));
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
W = sparse(1:meas.num.pmu, 1:meas.num.pmu, ...
           ones(meas.num.pmu, 1));
% -------------------------------------------------------------------------

% ------------------------- Construct H matrix --------------------------
[ rowIi, colIi ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));  

% ---------------------------- Row indices --------------------------------
iHIji = [ 1:meas.num.pIji, 1:meas.num.pIji ];
iHIij = [ (meas.num.pIji  + (1:meas.num.pIij)), (meas.num.pIji  + (1:meas.num.pIij)) ];
iHV = (meas.num.pIji + meas.num.pIij + (1:meas.num.pV));
iHIi = meas.num.pIji + meas.num.pIij + meas.num.pV + rowIi;
% -------------------------------------------------------------------------

% ------------------------ Column indices ---------------------------------
jHIji = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)); ];
jHIij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.Iij));
                 ];
jHV = meas.pmu.onbus(meas.pmu.v);
jHIi = colIi;
% -------------------------------------------------------------------------

% ------------------------- Values of elements ----------------------------
vHIji = [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)); ];
vHIij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); 
    ];
vHV = ones(meas.num.pV, 1);
vHIi = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Ii), :));
% -------------------------------------------------------------------------

H = sparse(...
             [ iHIji, iHIij, iHV, iHIi  ],... 
             [ jHIji; jHIij; jHV; jHIi  ], ...
             [ vHIji; vHIij; vHV; vHIi  ], ...
             meas.num.pmu, powsys.num.bus ...
         );
% -------------------------------------------------------------------------

% ---------------------- Construct J matrix -------------------------------
J = powsys.ybus.y(powsys.bus.zi, :);
% -------------------------------------------------------------------------

% --------------------------- Solve the LS problem ------------------------
x = [ H' * W * H,       J'
       J    sparse(powsys.num.zi, powsys.num.zi) ] \ ...
           [ H' *  W * z
                 sparse(powsys.num.zi, 1) ];
% -------------------------------------------------------------------------

% ----------------------- Estimator results -------------------------------
Vc = x;
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end