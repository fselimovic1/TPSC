function [ Vc, iter, converged, info ] = run_fEKFrect_rtse(powsys, meas, rtsesettings, Vc, initialStage)
% ------------ Matrices computed only at the inital run -------------------
if initialStage
    global H R Q;
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
    accI = accI + meas.num.pV; 
    iHvRe = (accI + (1:meas.num.pV));
    accI = accI + meas.num.pV; 
    iHvIm = (accI + (1:meas.num.pV));
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
    % ---------------------------------------------------------------------
    
    % ----------------------- Elements' values ----------------------------
    vHibrOre = [ real(powsys.ybus.tofrom(meas.pmu.loc(meas.pmu.ibranchO)));
                 -imag(powsys.ybus.tofrom(meas.pmu.loc(meas.pmu.ibranchO)));
                 real(powsys.ybus.toto(meas.pmu.loc(meas.pmu.ibranchO)));
                 -imag(powsys.ybus.toto(meas.pmu.loc(meas.pmu.ibranchO)))];
    vHibrOim = [ imag(powsys.ybus.tofrom(meas.pmu.loc(meas.pmu.ibranchO)));
                 real(powsys.ybus.tofrom(meas.pmu.loc(meas.pmu.ibranchO)));
                 imag(powsys.ybus.toto(meas.pmu.loc(meas.pmu.ibranchO)));
                 real(powsys.ybus.toto(meas.pmu.loc(meas.pmu.ibranchO)))];
    vHibrRe = [   real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 -imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)));
                 -imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)))];
    vHibrIm = [ imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)));
                 imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)));
                 real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)))];
    vHv = ones(meas.num.pV, 1);
    % ---------------------------------------------------------------------
    H = sparse([ iHibrOre, iHibrOim, iHibrRe, iHibrIm, iHvRe, iHvIm ],...
               [ jHibrO, jHibrO, jHibr, jHibr, jHvRe, jHvIm ],...
               [ vHibrOre, vHibrOim, vHibrRe, vHibrIm, vHv, vHv ]);
    R = 10^-9 .* eye(2 * meas.num.pmu);
    Q = 10^-9 .* eye(2 * powsys.num.bus);
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

% ------------------------ Measurements vector ----------------------------
z = [
      meas.pmu.m(meas.pmu.ibranchO) .* cos(meas.pmu.a(meas.pmu.ibranchO));
      meas.pmu.m(meas.pmu.ibranch) .* cos(meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* cos(meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.inj) .* cos(meas.pmu.a(meas.pmu.inj));
      meas.pmu.m(meas.pmu.ibranchO) .* sin(meas.pmu.a(meas.pmu.ibranchO));
      meas.pmu.m(meas.pmu.ibranch) .* sin(meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* sin(meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.inj) .* sin(meas.pmu.a(meas.pmu.inj));
    ];
% -------------------------------------------------------------------------
end

