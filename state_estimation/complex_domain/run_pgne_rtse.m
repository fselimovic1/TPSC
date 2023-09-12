function [ Vc, iter, converged, info ] = run_pgne_rtse(powsys, meas, rtsesettings, Vc)
info.method = 'Complex Perturbed Gauss-Newton Estimator';
info.paper = 'A Complex Variable Perturbed Gauss-Newton Method for Tracking Mode State Estimation';
% -------------------- Eq. for the slack bus (no PMU case) ----------------
nopmu =  ~meas.num.pmu;
% -------------------------------------------------------------------------

% --------------------- Setting additional variables ----------------------
global L D H W zPh zPhKeys zS zSKeys;
iter = 1;
converged = 0;
% -------------------------------------------------------------------------

% --------------------- Initialize state variables ------------------------
x = Vc;
% -------------------------------------------------------------------------

% ------------ Matrices computed only at the inital run -------------------
if rtsesettings.initialStage
    zPh = zeros(meas.num.pmu, 1);
    zPhKeys = [ -meas.pmu.loc(meas.pmu.IijO)
                 powsys.num.branch  + meas.pmu.loc(meas.pmu.Iij)
                 2 * powsys.num.branch +  meas.pmu.loc(meas.pmu.v)
                 2 * powsys.num.branch + powsys.num.bus + meas.pmu.loc(meas.pmu.inj);
            ];
    % ---------------------------------------------------------------------

    % ------------------------ SCADA measurement keys ---------------------
    [ idxPijO, ~ ] = myintersect(meas.scada.loc, meas.scada.pijO, meas.scada.qijO);
    [ idxPij, ~ ] = myintersect(meas.scada.loc, meas.scada.pij, meas.scada.qij);
    [ idxPinj, ~ ] = myintersect(meas.scada.loc, meas.scada.pinj, meas.scada.qinj); 
    zS = zeors(numel(locSijO) + locSij + locSinj + meas.num.sIijmO + meas.num.sIijm + meas.num.sVm, 1);
    zSKeys = [ 
                -meas.scada.loc(idxPijO); meas.scada.loc(idxPij) + powsys.num.branch; ...
                meas.scada.loc(idxPinj) + 2 * powsys.num.branch;...
                -meas.scada.loc(meas.scada.IijmO) + 2 * powsys.num.branch + powsys.num.bus; ...
                meas.scada.loc(meas.scada.Iijm) + 3 * powsys.num.branch + powsys.num.bus; ...
                meas.scada.loc(meas.scada.vm) + 4 * powsys.num.branch + powsys.num.bus; ...
        ];  
    % ---------------------------------------------------------------------
    % --------------- Measurements Coefficient matrix - H -----------------
    [ rowInj, colInj ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));  
    % ---------------------------- Row indices ----------------------------
    iH11IijO = [ 1:meas.num.pIijO, 1:meas.num.pIijO ];
    iH11Iij = [ (meas.num.pIijO  + (1:meas.num.pIij)), (meas.num.pIijO  + (1:meas.num.pIij)) ];
    iH11V = (meas.num.pIijO + meas.num.pIij + (1:meas.num.pV));
    iH11Inj = meas.num.pIijO + meas.num.pIij + meas.num.pV + rowInj;
    % ---------------------------------------------------------------------

    % ------------------------ Column indices -----------------------------
    jH11Iij0 = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.IijO)); 
                 powsys.branch.j(-meas.pmu.loH(meas.pmu.IijO)); ];
    jH11Iij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.Iij));
                 ];
    jH11V = meas.pmu.onbus(meas.pmu.v);
    jH11Inj = colInj;
    % ---------------------------------------------------------------------

    % ------------------------- Values of elements ------------------------
    vH11IijO = [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.IijO));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.IijO)); ];
    vH11Iij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); 
    ];
    vH11V = ones(meas.num.pV, 1);
    vH11Inj = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));
    % ---------------------------------------------------------------------

    H = sparse(...
             [ iH11IijO, iH11Iij, iH11V, iH11Inj  ],... 
             [ jH11Iij0; jH11Iij; jH11V; jH11Inj  ], ...
             [ vH11IijO; vH11Iij; vH11V; vH11Inj  ], ...
             meas.num.pmu, powsys.num.bus...
         );
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------


% ------------------------ Adjust measurement set -------------------------
% -------------------------- PMU data -------------------------------------
[ ~, idxIijO ] = find(-meas.pmu.loc(meas.pmu.IijO), zPhKeys);
[ ~, idxIij ] = find(meas.pmu.loc(meas.pmu.IijO) + powsys.num.branch, zPhKeys);
[ ~, idxV ] = find(meas.pmu.loc(meas.pmu.v)  + 2 * powsys.num.branch, zPhKeys);
[ ~, idxIinj ] = find(meas.pmu.loc(meas.pmu.Iinj) + 2 * powsys.num.branch + powsys.num.bus, zPhKeys);
zPh(idxIijO, idxIij, idxV, idxIinj) = [ meas.pmu.m(meas.pmu.IijO) .* exp(1i .* meas.pmu.a(meas.pmu.IijO));
        meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
        meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
        meas.pmu.m(meas.pmu.Iinj) .* exp(1i .* meas.pmu.a(meas.pmu.Iinj)); ];
% -------------------------------------------------------------------------

% ---------------------------- SCADA data ---------------------------------
[ idxPijO, idxQijO ] = myintersect(meas.scada.loc, meas.scada.pijO, meas.scada.qijO);
[ idxPij, idxQij ] = myintersect(meas.scada.loc, meas.scada.pij, meas.scada.qij);
[ idxPinj, idxQinj ] = myintersect(meas.scada.loc, meas.scada.pinj, meas.scada.qinj); 
[ ~, idxIijmO ] = find(meas.pmu.loc(meas.pmu.IijmO) + 2 * powsys.num.branch + powsys.num.bus, zSKeys);
[ ~, idxIijm ] = find(meas.pmu.loc(meas.pmu.Iijm)  + 3 * powsys.num.branch + powsys.num.bus, zSKeys);
[ ~, idxVm ] = find(meas.pmu.loc(meas.pmu.vm) + 4 * powsys.num.branch + powsys.num.bus, zSKeys);
zS([-idxPijO; idxPij + powsys.num.branch; idxPinj + 2 * powsys.num.branch;...
    idxIijmO + 2 * powsys.num.branch + powsys.num.bus;  ...
    idxIijm + 3 * powsys.num.branch + powsys.num.bus; ...
    idxVm + 4 * powsys.num.branch + powsys.num.bus]) = [
    conj(meas.scada.m(idxPijO) + 1i .* meas.scada.m(idxQijO)) ...
    ./ conj(x(meas.scada.bus(idxPijO)));
    conj(meas.scada.m(idxPij) + 1i .* meas.scada.m(idxQij)) ...
    ./ conj(x(meas.scada.bus(idxPij)));
    conj(meas.scada.m(idxPinj) + 1i .* meas.scada.m(idxQinj)) ...
    ./ conj(meas.scada.bus(idxPinj));    
    meas.scada.m(meas.scada.IijmO) .* exp(1i .* 5);
    meas.scada.m(meas.scada.Iijm);
    meas.scada.m(meas.scada.vm)
  ];
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
end

