function [ Vc, iter, converged, info ] = run_pgne_rtse(powsys, meas, rtsesettings, Vc)
info.method = 'Complex Perturbed Gauss-Newton Estimator';
info.paper = 'A Complex Variable Perturbed Gauss-Newton Method for Tracking Mode State Estimation';
% -------------------- Eq. for the slack bus (no PMU case) ----------------
nopmu =  ~meas.num.pmu;
% -------------------------------------------------------------------------

% --------------------- Setting additional variables ----------------------
global L D H W;
iter = 1;
converged = 0;
% -------------------------------------------------------------------------

% --------------------- Initialize state variables ------------------------
if nargin == 4
    x = varargin{1};
else
    if sesettings.flatStart
        x = ones( 2 * powsys.num.bus, 1);
    else
        x = [ powsys.bus.Vmi .* exp(1i * powsys.bus.Vai);
              powsys.bus.Vmi .* exp(-1i * powsys.bus.Vai);
             ];
    end
end
% -------------------------------------------------------------------------

% ------------ Matrices computed only at the inital run -------------------
if rtsesettings.initialStage
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

end

