function [ x, x_, converged, info ] = run_dkfQA_dse(powsys, meas, dsesettings, x_)
info.method = "Discrete Kalman Filter with the online assessment of Q";
info.paper = [ 'Performance Assessment of Linear State Estimators ', ...
               'Using Synchrophasor Measurements'...
               ];
% ------------------------- Variables definition --------------------------
converged = 1;
global H R P_;
reidx = 1:2:2 * powsys.num.bus - 1;
imidx = 2:2:2 * powsys.num.bus;
Trr = 1 / dsesettings.fc;
% -------------------------------------------------------------------------

% ------------ Matrices computed only at the inital run -------------------
if dsesettings.initialStage
    % -------------------- Measurement matrix - H -------------------------
    accI = 0;
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
    % ---------------------------------------------------------------------
    H = sparse([ iHIjiRe, iHIjiIm, iHIijRe, iHIijIm, iHvRe, iHvIm, ],...
               [ jHIji; jHIji; jHIij; jHIij; jHvRe; jHvIm; ],...
               [ vHIjiRe; vHIjiIm; vHIijRe; vHIijIm; vHv; vHv; ]);
    % ---------------------------------------------------------------------
    
    
    % ------------------ Measurement noise covarinace matrix --------------
    
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------

end

