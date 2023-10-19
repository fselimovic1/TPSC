function [ Vc, iter, converged, info ] = run_wls_rect_sse(powsys, meas)
info.method = "WLS with state variables in rectangular coordinates";
global weight;
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
iRIji = [ (1:meas.num.pIji), (1:meas.num.pIji),... 
           accI + (1:meas.num.pIji), accI + (1:meas.num.pIji) ];    
accI = accI + meas.num.pIji;
iRIij = [ accI + (1:meas.num.pIij), accI + (1:meas.num.pIij),... 
          accI + meas.num.pIij + (1:meas.num.pIij), accI + meas.num.pIij + (1:meas.num.pIij) ];
accI = accI + 2 * meas.num.pIij;
iRv = [ accI + (1:meas.num.pV), accI + (1:meas.num.pV), ...
        accI + meas.num.pV + (1:meas.num.pV), accI + meas.num.pV + (1:meas.num.pV) ];
iRidx = [ iRIji, iRIij, iRv ];
% -------------------------------------------------------------------------
    
% ---------------------------- Column indices -----------------------------
accJ = meas.num.pIji;
jRIji = [ 1:meas.num.pIji, (1:meas.num.pIji) + meas.num.pIji, ...
           1:meas.num.pIji, accJ + (1:meas.num.pIji) ];
accJ = accJ + meas.num.pIji;
jRIij = [ (accJ + (1:meas.num.pIij)), (accJ + (1:meas.num.pIij)) + meas.num.pIij,... 
           (accJ + (1:meas.num.pIij)), (accJ + (1:meas.num.pIij)) + meas.num.pIij ];
accJ = accJ + 2 * meas.num.pIij;   
jRv = [ (accJ + (1:meas.num.pV)), (accJ + (1:meas.num.pV)) + meas.num.pV,... 
         (accJ + (1:meas.num.pV)) + meas.num.pV, (accJ + (1:meas.num.pV)) ];
jRidx = [ jRIji, jRIij, jRv ];
% -------------------------------------------------------------------------
% ----------------- Compute values of the elements ------------------------
vRIji = [ (cos(meas.pmu.a(meas.pmu.Iji)).^2) .* (meas.pmu.msd(meas.pmu.Iji) .^2)...
            + meas.pmu.m(meas.pmu.Iji).^2 .* (sin(meas.pmu.a(meas.pmu.Iji)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);
            cos(meas.pmu.a(meas.pmu.Iji)) .* sin(meas.pmu.a(meas.pmu.Iji)) .* ...
            meas.pmu.msd(meas.pmu.Iji) .^2  - meas.pmu.m(meas.pmu.Iji).^2 .* ...
            cos(meas.pmu.a(meas.pmu.Iji)) .* sin(meas.pmu.a(meas.pmu.Iji)) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);
            cos(meas.pmu.a(meas.pmu.Iji)) .* sin(meas.pmu.a(meas.pmu.Iji)) .* ...
            meas.pmu.msd(meas.pmu.Iji) .^2  - meas.pmu.m(meas.pmu.Iji).^2 .* ...
            cos(meas.pmu.a(meas.pmu.Iji)) .* sin(meas.pmu.a(meas.pmu.Iji)) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);
            (sin(meas.pmu.a(meas.pmu.Iji)).^2) .* (meas.pmu.msd(meas.pmu.Iji) .^2)...
            + meas.pmu.m(meas.pmu.Iji).^2 .* (cos(meas.pmu.a(meas.pmu.Iji)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iji) .^2);
             ];      
vRIij = [ (cos(meas.pmu.a(meas.pmu.Iij)).^2) .* (meas.pmu.msd(meas.pmu.Iij) .^2)...
            + meas.pmu.m(meas.pmu.Iij).^2 .* (sin(meas.pmu.a(meas.pmu.Iij)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);
            cos(meas.pmu.a(meas.pmu.Iij)) .* sin(meas.pmu.a(meas.pmu.Iij)) .* ...
            meas.pmu.msd(meas.pmu.Iij) .^2  - meas.pmu.m(meas.pmu.Iij).^2 .* ...
            cos(meas.pmu.a(meas.pmu.Iij)) .* sin(meas.pmu.a(meas.pmu.Iij)) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);
            cos(meas.pmu.a(meas.pmu.Iij)) .* sin(meas.pmu.a(meas.pmu.Iij)) .* ...
            meas.pmu.msd(meas.pmu.Iij) .^2  - meas.pmu.m(meas.pmu.Iij).^2 .* ...
            cos(meas.pmu.a(meas.pmu.Iij)) .* sin(meas.pmu.a(meas.pmu.Iij)) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);
            (sin(meas.pmu.a(meas.pmu.Iij)).^2) .* (meas.pmu.msd(meas.pmu.Iij) .^2)...
            + meas.pmu.m(meas.pmu.Iij).^2 .* (cos(meas.pmu.a(meas.pmu.Iij)).^2) .* ...
            (meas.pmu.asd(meas.pmu.Iij) .^2);
             ]; 
vRv = [ (cos(meas.pmu.a(meas.pmu.v)).^2) .* (meas.pmu.msd(meas.pmu.v) .^2)...
            + meas.pmu.m(meas.pmu.v).^2 .* (sin(meas.pmu.a(meas.pmu.v)).^2) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);
            cos(meas.pmu.a(meas.pmu.v)) .* sin(meas.pmu.a(meas.pmu.v)) .* ...
            meas.pmu.msd(meas.pmu.v) .^2  - meas.pmu.m(meas.pmu.v).^2 .* ...
            cos(meas.pmu.a(meas.pmu.v)) .* sin(meas.pmu.a(meas.pmu.v)) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);
            cos(meas.pmu.a(meas.pmu.v)) .* sin(meas.pmu.a(meas.pmu.v)) .* ...
            meas.pmu.msd(meas.pmu.v) .^2  - meas.pmu.m(meas.pmu.v).^2 .* ...
            cos(meas.pmu.a(meas.pmu.v)) .* sin(meas.pmu.a(meas.pmu.v)) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);
            (sin(meas.pmu.a(meas.pmu.v)).^2) .* (meas.pmu.msd(meas.pmu.v) .^2)...
            + meas.pmu.m(meas.pmu.v).^2 .* (cos(meas.pmu.a(meas.pmu.v)).^2) .* ...
            (meas.pmu.asd(meas.pmu.v) .^2);
             ];    
vR = [  1 ./ vRIji; 1 ./ vRIij; 1 ./ vRv; ];
% -------------------------------------------------------------------------
R = sparse(iRidx, jRidx, vR);
% -------------------------------------------------------------------------

% ---------------------- Measurement matrix - H ---------------------------
accI = 0;
% -------------------------- Row indices ----------------------------------
iHIjire = [ 1:meas.num.pIji, 1:meas.num.pIji, 1:meas.num.pIji, 1:meas.num.pIji ];
accI = accI + meas.num.pIji;
iHIjiim = [ (accI + (1:meas.num.pIji)), (accI + (1:meas.num.pIji)),...
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
% -------------------------------------------------------------------------

% ---------------------------- Column indices -----------------------------
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
% -------------------------------------------------------------------------

% ------------------------- Elements' values ------------------------------
vHIjire = [ real(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
             -imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
             real(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)));
             -imag(powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)))];
vHIjiim = [ imag(powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji)));
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
% -------------------------------------------------------------------------
H = sparse([ iHIjire, iHIjiim, iHIijRe, iHIijIm, iHvRe, iHvIm ],...
           [ jHIji; jHIji; jHIij; jHIij; jHvRe; jHvIm ],...
           [ vHIjire; vHIjiim; vHIijRe; vHIijIm; vHv; vHv ]);
% -------------------------------------------------------------------------

% --------------------------- Solve the LS problem ------------------------
if ~weight
    R = sparse(1:2 * meas.num.pmu, 1:2 * meas.num.pmu, 1);
end
x = (H' * R * H) \ (H' * R * z);
% -------------------------------------------------------------------------

% ----------------------- Estimator results -------------------------------
Vc = x;
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end

