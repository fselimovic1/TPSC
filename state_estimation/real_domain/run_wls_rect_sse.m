function [ Vc, iter, converged, info ] = run_wls_rect_sse(powsys, meas)
info.method = "WLS with state variables in rectangular coordinates";

% --------------------- Setting additional variables ----------------------
iter = 1;
converged = 1;
% -------------------------------------------------------------------------


% ------------------------ Measurements vector ----------------------------
z = [
      meas.pmu.m(meas.pmu.Iji) .* cos(meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iji) .* sin(meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* cos(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.Iij) .* sin(meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* cos(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.v) .* sin(meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* cos(meas.pmu.a(meas.pmu.Iinj));
      meas.pmu.m(meas.pmu.Iinj) .* sin(meas.pmu.a(meas.pmu.Iinj));
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
% ----------------------------- Row indices -------------------------------
accI = meas.num.pIji;
iRIbrO = [ (1:meas.num.pIji),... 
           accI + (1:meas.num.pIji)];    
accI = accI + meas.num.pIji;
iRIbr = [ accI + (1:meas.num.pIij),... 
          accI + meas.num.pIij + (1:meas.num.pIij)];
accI = accI + 2 * meas.num.pIij;
iRv = [ accI + (1:meas.num.pV),...
        accI + meas.num.pV + (1:meas.num.pV)];
iRidx = [ iRIbrO, iRIbr, iRv ];
% -------------------------------------------------------------------------
    
% ---------------------------- Column indices -----------------------------
accJ = meas.num.pIji;
jRIbrO = [ 1:meas.num.pIji,...
           accJ + (1:meas.num.pIji) ];
accJ = accJ + meas.num.pIji;
jRIbr = [ (accJ + (1:meas.num.pIij)),... 
           accJ + meas.num.pIij + (1:meas.num.pIij)];
accJ = accJ + 2 * meas.num.pIij;   
jRv = [ (accJ + (1:meas.num.pV)),... 
         accJ + meas.num.pV + (1:meas.num.pV)];
jRidx = [ jRIbrO, jRIbr, jRv ];
% -------------------------------------------------------------------------
% ----------------- Compute values of the elements ------------------------
vRIbrO = [ (cos(meas.pmu.a(meas.pmu.Iji)).^2) .* (meas.pmu.msd(meas.pmu.Iji) .^2)...
            + meas.pmu.m(meas.pmu.Iji).^2 .* (sin(meas.pmu.a(meas.pmu.Iji)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);...
            (sin(meas.pmu.a(meas.pmu.Iji)).^2) .* (meas.pmu.msd(meas.pmu.Iji) .^2)...
            + meas.pmu.m(meas.pmu.Iji).^2 .* (cos(meas.pmu.a(meas.pmu.Iji)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);
             ];      
vRIbr = [ (cos(meas.pmu.a(meas.pmu.Iij)).^2) .* (meas.pmu.msd(meas.pmu.Iij) .^2)...
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
vR = [  1 ./ vRIbrO; 1 ./ vRIbr; 1 ./ vRv; ];
% -------------------------------------------------------------------------
R = sparse(iRidx, jRidx, vR);
% -------------------------------------------------------------------------

% ---------------------- Measurement matrix - H ---------------------------
accI = 0;
% -------------------------- Row indices ----------------------------------
iHibrOre = [ 1:meas.num.pIji, 1:meas.num.pIji, 1:meas.num.pIji, 1:meas.num.pIji ];
accI = accI + meas.num.pIji;
iHibrOim = [ (accI + (1:meas.num.pIji)), (accI + (1:meas.num.pIji)),...
            (accI + (1:meas.num.pIji)), (accI + (1:meas.num.pIji)) ];
accI = accI + meas.num.pIji; 
iHibrRe = [ (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)),...
            (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)) ];
accI = accI + meas.num.pIij;         
iHibrIm = [ (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)),...
            (accI + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)) ]; 
accI = accI + meas.num.pIij; 
iHvRe = (accI + (1:meas.num.pV));
accI = accI + meas.num.pV; 
iHvIm = (accI + (1:meas.num.pV));
% -------------------------------------------------------------------------

% ---------------------------- Column indices -----------------------------
jHibrO = [ 2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji)) - 1; 
           2 .* powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji));
           2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)) - 1;
           2 .* powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)) ];
jHibr = [  2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)) - 1; 
           2 .* powsys.branch.i(meas.pmu.loc(meas.pmu.Iij));
           2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.Iij)) - 1;
           2 .* powsys.branch.j(meas.pmu.loc(meas.pmu.Iij)) ];
jHvRe  = 2 .* meas.pmu.onbus(meas.pmu.v) - 1;
jHvIm  = 2 .* meas.pmu.onbus(meas.pmu.v);
% -------------------------------------------------------------------------

% ------------------------- Elements' values ------------------------------
vHibrOre = [ real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
             -imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
             real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)));
             -imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)))];
vHibrOim = [ imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
             real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
             imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)));
             real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)))];
vHibrRe = [   real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
             -imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
             real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
             -imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
vHibrIm = [ imag(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
             real(powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij)));
             imag(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)));
             real(powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)))];
vHv = ones(meas.num.pV, 1);
% -------------------------------------------------------------------------
H = sparse([ iHibrOre, iHibrOim, iHibrRe, iHibrIm, iHvRe, iHvIm ],...
           [ jHibrO; jHibrO; jHibr; jHibr; jHvRe; jHvIm ],...
           [ vHibrOre; vHibrOim; vHibrRe; vHibrIm; vHv; vHv ]);
% -------------------------------------------------------------------------

% --------------------------- Solve the LS problem ------------------------
% R = sparse(1:2 * meas.num.pmu, 1:2 * meas.num.pmu, 1);
x = (H' * R * H) \ (H' * R * z);
% -------------------------------------------------------------------------

% ----------------------- Estimator results -------------------------------
Vc = x;
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end

