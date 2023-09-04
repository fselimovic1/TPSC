function [ Vc, iter, converged, info ] = run_cls_sse(powsys, meas)
info.method = 'Linear State Estimation (PMU only) in Complex Variables';
% --------------------- Setting additional variables ----------------------
iter = 1;
converged = 1;
% -------------------------------------------------------------------------

% ----------------------- Vector of measurement values --------------------
z = [ meas.pmu.m(meas.pmu.IijO) .* exp(1i .* meas.pmu.a(meas.pmu.IijO));
      meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* exp(1i .* meas.pmu.a(meas.pmu.Iinj));
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
% W = sparse(1:meas.num.pmu, 1:meas.num.pmu, ...
%            [ str2double(sesettings.mweights(2)) .* ones( 2 * meas.num.pmu, 1);...
%            ones(meas.num.scada, 1) ]);
% -------------------------------------------------------------------------

% ------------------------- Construct H matrix --------------------------
[ rowInj, colInj ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));  

% ---------------------------- Row indices --------------------------------
iHIijO = [ 1:meas.num.pIijO, 1:meas.num.pIijO ];
iHIij = [ (meas.num.pIijO  + (1:meas.num.pIij)), (meas.num.pIijO  + (1:meas.num.pIij)) ];
iHV = (meas.num.pIijO + meas.num.pIij + (1:meas.num.pV));
iHInj = meas.num.pIijO + meas.num.pIij + meas.num.pV + rowInj;
% -------------------------------------------------------------------------

% ------------------------ Column indices ---------------------------------
jHIij0 = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.IijO)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.IijO)); ];
jHIij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.Iij));
                 ];
jHV = meas.pmu.onbus(meas.pmu.v);
jHInj = colInj;
% -------------------------------------------------------------------------

% ------------------------- Values of elements ----------------------------
vHIijO = [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.IijO));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.IijO)); ];
vHIij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); 
    ];
vHV = ones(meas.num.pV, 1);
vHInj = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));
% -------------------------------------------------------------------------

H = sparse(...
             [ iHIijO, iHIij, iHV, iHInj  ],... 
             [ jHIij0; jHIij; jHV; jHInj  ], ...
             [ vHIijO; vHIij; vHV; vHInj  ], ...
             meas.num.pmu, powsys.num.bus ...
         );
% -------------------------------------------------------------------------


% --------------------------- Solve the LS problem ------------------------
    x = H \ z;
% -------------------------------------------------------------------------

% ----------------------- Estimator results -------------------------------
Vc = x;
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end