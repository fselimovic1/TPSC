function [ Vc, iter, converged, info ] = run_pgne_tse(powsys, meas, tsesettings, Vc)
info.method = 'Complex Perturbed Gauss-Newton Estimator';
info.paper = 'A Complex Variable Perturbed Gauss-Newton Method for Tracking Mode State Estimation';
info.nonZerosInH = 0;
info.redundancy = 0;
% -------------------- Eq. for the slack bus (no PMU case) ----------------
nopmu =  ~meas.num.pmu;
% -------------------------------------------------------------------------

% --------------------- Setting additional variables ----------------------
global L D H W zPh zPhKeys zS zSKeys nEQ;
iter = 1;
converged = 0;
% -------------------------------------------------------------------------

% --------------------- Initialize state variables ------------------------
x = Vc;
% -------------------------------------------------------------------------

% ------------ Matrices computed only at the inital run -------------------
    %---------------------------------------------------------------------
if tsesettings.initialStage
    %------------------------ PMU measurement keys ------------------------
    zPh = zeros(meas.num.pmu, 1);
    zPhKeys = [ -meas.pmu.loc(meas.pmu.IijO)
                 powsys.num.branch  + meas.pmu.loc(meas.pmu.Iij)
                 2 * powsys.num.branch +  meas.pmu.loc(meas.pmu.v)
                 2 * powsys.num.branch + powsys.num.bus + meas.pmu.loc(meas.pmu.Iinj);
            ];
    % ---------------------------------------------------------------------

    % ------------------------ SCADA measurement keys ---------------------
    [ idxPijO, ~ ] = myintersect(meas.scada.loc, meas.scada.pijO, meas.scada.qijO);
    [ idxPij, ~ ] = myintersect(meas.scada.loc, meas.scada.pij, meas.scada.qij);
    [ idxPinj, ~ ] = myintersect(meas.scada.loc, meas.scada.pinj, meas.scada.qinj); 
    zS = zeros(numel(idxPijO) + numel(idxPij) + numel(idxPinj) +...
         meas.num.sIijmO + meas.num.sIijm + meas.num.sVm, 1);
    zSKeys = [ 
                -meas.scada.loc(idxPijO); meas.scada.loc(idxPij) + powsys.num.branch; ...
                meas.scada.loc(idxPinj) + 2 * powsys.num.branch;...
                -meas.scada.loc(meas.scada.IijmO) + 2 * powsys.num.branch + powsys.num.bus; ...
                meas.scada.loc(meas.scada.Iijm) + 3 * powsys.num.branch + powsys.num.bus; ...
                meas.scada.loc(meas.scada.vm) + 4 * powsys.num.branch + powsys.num.bus; ...
        ];  
    % ---------------------------------------------------------------------
    
    
    % --------------- Measurements Coefficient matrix - H -----------------
    % ----------------------- PMU measurements ----------------------------
    [ rowInj, colInj ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));  
    % ---------------------------- Row indices ----------------------------
    accI = 0;
    iH11IijO = [ 1:meas.num.pIijO, 1:meas.num.pIijO ];
    accI = accI + meas.num.pIijO;
    iH11Iij = [ (accI  + (1:meas.num.pIij)), (accI + (1:meas.num.pIij)) ];
    accI = accI + meas.num.pIij;
    iH11V = (accI + (1:meas.num.pV));
    accI = accI + meas.num.pV;
    iH11Inj = accI + rowInj;
    if ~isempty(rowInj)
        accI = accI + max(rowInj);
    end
    % ---------------------------------------------------------------------

    % ------------------------ Column indices -----------------------------
    jH11IijO = [    powsys.branch.i(-meas.pmu.loc(meas.pmu.IijO)); 
                    powsys.branch.j(-meas.pmu.loc(meas.pmu.IijO)); 
                    ];
    jH11Iij =  [    powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
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
     
     % ---------------------- Adapted SCADA measurements ------------------
    [ rowInj, colInj ] = find(powsys.ybus.y(meas.scada.loc(idxPinj), :));  
    % ---------------------------- Row indices ----------------------------
    iH21SijO = [ accI + (1:numel(idxPijO)), accI + (1:numel(idxPijO)) ];
    accI = accI + numel(idxPijO);
    iH21Sij = [ (accI  + (1:numel(idxPij))), (accI  + (1:numel(idxPij))) ];
    accI = accI + numel(idxPij);
    iH21Sinj = accI + rowInj;
    if ~isempty(rowInj)
        accI = accI + max(rowInj);
    end
    iH21IijmO = accI + [ (1:meas.num.sIijmO), (1:meas.num.sIijmO) ];
    accI = accI + meas.num.sIijmO;
    iH21Iijm = accI + [ (1:meas.num.sIijm), (1:meas.num.sIijm) ];
    accI = accI + meas.num.sIijm;
    iH21Vm = accI + (1:meas.num.sVm);
    accI = accI + meas.num.sVm;
    numScada = accI - meas.num.pmu;
    % ---------------------------------------------------------------------

    % ------------------------ Column indices -----------------------------
    jH21SijO = [ 
                 powsys.branch.i(-meas.scada.loc(idxPijO)); 
                 powsys.branch.j(-meas.scada.loc(idxPijO));
                 ];
    jH21Sij = [ 
                 powsys.branch.i(meas.scada.loc(idxPij)); 
                 powsys.branch.j(meas.scada.loc(idxPij));
                 ];
    jH21Sinj = colInj;
    jH21IijmO = [ 
                 powsys.branch.i(-meas.scada.loc(meas.scada.IijmO)); 
                 powsys.branch.j(-meas.scada.loc(meas.scada.IijmO));
                 ];
    jH21Iijm = [ 
                 powsys.branch.i(meas.scada.loc(meas.scada.Iijm)); 
                 powsys.branch.j(meas.scada.loc(meas.scada.Iijm));
                 ];
    jH21Vm = meas.scada.onbus(meas.scada.vm);
    % ---------------------------------------------------------------------

    % ------------------------- Values of elements ------------------------
    vH21SijO = [ 
                 powsys.ybus.tofrom(-meas.scada.loc(idxPijO));...
                 powsys.ybus.toto(-meas.scada.loc(idxPijO));
                 ];
    vH21Sij = [ 
                 powsys.ybus.fromfrom(meas.scada.loc(idxPij));...
                 powsys.ybus.fromto(meas.scada.loc(idxPij));
                 ];
    vH21Sinj = nonzeros(powsys.ybus.y(meas.scada.loc(idxPinj), :));
    vH21IijmO = [ 
                 powsys.ybus.tofrom(-meas.scada.loc(meas.scada.IijmO));...
                 powsys.ybus.toto(-meas.scada.loc(meas.scada.IijmO));
                 ];
    vH21Iijm = [ 
                 powsys.ybus.fromfrom(meas.scada.loc(meas.scada.Iijm));...
                 powsys.ybus.fromto(meas.scada.loc(meas.scada.Iijm));
                 ];
    vH21Vm = ones(meas.num.sVm, 1);
    % ---------------------------------------------------------------------
    H = sparse(...
             [ iH11IijO, iH11Iij, iH11V, iH11Inj, iH21SijO, iH21Sij, iH21Sinj', ...
                iH21IijmO, iH21Iijm, iH21Vm ],... 
             [ jH11IijO; jH11Iij; jH11V; jH11Inj; jH21SijO; jH21Sij; jH21Sinj; ...
                jH21IijmO; jH21Iijm; jH21Vm  ], ...
             [ vH11IijO; vH11Iij; vH11V; vH11Inj; vH21SijO; vH21Sij; vH21Sinj; ...
                vH21IijmO; vH21Iijm; vH21Vm  ] );
    % ---------------------------------------------------------------------
    
    % ---------------------- Measurement weights --------------------------
    wIdx = 1:(meas.num.pmu + numScada);
    W = sparse(wIdx, wIdx, [ str2double(tsesettings.mweights(2)) .* ones( meas.num.pmu, 1);...
           ones(numScada, 1) ]);
    % ---------------------------------------------------------------------
    % ---------------------- Construct J matrix ---------------------------
    rJzi = powsys.ybus.y(powsys.bus.zi, :);
    if nopmu
        nEQ = powsys.num.zi + 1;
        rJsl = sparse(1, powsys.num.islack, -1i / 2, 1, powsys.num.bus);
        J = [rJzi; rJsl];
    else
        nEQ = powsys.num.zi;
        J =  rJzi;
    end
    % ---------------------------------------------------------------------
    
    % ----------------- Factorize system matrix ---------------------------
    A = full([   H' * W * H        J';
                        J      sparse(nEQ, nEQ); ]);
    A = A - 1i .* imag(diag(diag(A)));
    [ L, D ] = ldl(A);
    L = sparse(L);
    D = sparse(D);
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------


% ------------------------ Adjust measurement set -------------------------
% -------------------------- PMU data -------------------------------------
[ ~, idxIijO ] = ismember(-meas.pmu.loc(meas.pmu.IijO), zPhKeys);
[ ~, idxIij ] = ismember(meas.pmu.loc(meas.pmu.Iij) + powsys.num.branch, zPhKeys);
[ ~, idxV ] = ismember(meas.pmu.loc(meas.pmu.v)  + 2 * powsys.num.branch, zPhKeys);
[ ~, idxIinj ] = ismember(meas.pmu.loc(meas.pmu.Iinj) + 2 * powsys.num.branch + powsys.num.bus, zPhKeys);
zPh([idxIijO; idxIij; idxV; idxIinj]) = [ ...
        meas.pmu.m(meas.pmu.IijO) .* exp(1i .* meas.pmu.a(meas.pmu.IijO));
        meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
        meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
        meas.pmu.m(meas.pmu.Iinj) .* exp(1i .* meas.pmu.a(meas.pmu.Iinj)); 
        ];
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

while iter < tsesettings.maxNumberOfIter 
    % ---------------------------- SCADA data ---------------------------------
    [ idxPijO, idxQijO ] = myintersect(meas.scada.loc, meas.scada.pijO, meas.scada.qijO);
    [ ~, idxSijO ] = ismember(-meas.scada.loc(idxPijO), zSKeys);
    [ idxPij, idxQij ] = myintersect(meas.scada.loc, meas.scada.pij, meas.scada.qij);
    [ ~, idxSij ] = ismember(meas.scada.loc(idxPij) + powsys.num.branch, zSKeys);
    [ idxPinj, idxQinj ] = myintersect(meas.scada.loc, meas.scada.pinj, meas.scada.qinj);
    [ ~, idxSinj ]  = ismember(meas.scada.loc(idxPinj) + 2 * powsys.num.branch, zSKeys) ;
    [ ~, idxIijmO ] = ismember(-meas.scada.loc(meas.scada.IijmO) + 2 * powsys.num.branch + powsys.num.bus, zSKeys);
    [ ~, idxIijm ] = ismember(meas.scada.loc(meas.scada.Iijm)  + 3 * powsys.num.branch + powsys.num.bus, zSKeys);
    [ ~, idxVm ] = ismember(meas.scada.loc(meas.scada.vm) + 4 * powsys.num.branch + powsys.num.bus, zSKeys);

    % ------------------------ Compute branch currents  -------------------
    locIijm = meas.scada.loc(meas.scada.Iijm);
    Iijm = sum([ powsys.ybus.fromfrom(locIijm), powsys.ybus.fromto(locIijm) ]...
                 .* [ Vc(powsys.branch.i(locIijm)), Vc(powsys.branch.j(locIijm))], 2);
    locIijmO = -meas.scada.loc(meas.scada.IijmO);
    IijmO = sum([ powsys.ybus.tofrom(locIijmO), powsys.ybus.toto(locIijmO) ]...
                .* [ Vc(powsys.branch.i(locIijmO)), Vc(powsys.branch.j(locIijmO))], 2);
    % ---------------------------------------------------------------------

    zS([idxSijO; idxSij; idxSinj; idxIijmO; idxIijm; idxVm]) = [
            conj(meas.scada.m(idxPijO) + 1i .* meas.scada.m(idxQijO)) ...
            ./ conj(x(meas.scada.onbus(idxPijO)));
            conj(meas.scada.m(idxPij) + 1i .* meas.scada.m(idxQij)) ...
            ./ conj(x(meas.scada.onbus(idxPij)));
            conj(meas.scada.m(idxPinj) + 1i .* meas.scada.m(idxQinj)) ...
            ./ conj(x(meas.scada.onbus(idxPinj)));    
            meas.scada.m(meas.scada.IijmO) .* exp(1i .* angle(IijmO));
            meas.scada.m(meas.scada.Iijm) .* exp(1i .* angle(Iijm));
            meas.scada.m(meas.scada.vm) .* exp(1i .* angle(x(meas.scada.onbus(meas.scada.vm))));
      ];
    % ---------------------------------------------------------------------
    xnew = L' \ ( D \ ( L \ [...
           H' * W * [ zPh; zS ]
           sparse(nEQ, 1)
                ]));
    % ---------------------------------------------------------------------         
    % ------------------------ Check Convergence --------------------------
    if max(abs(xnew(powsys.bus.busnew) - x)) < tsesettings.eps
        x = xnew(powsys.bus.busnew);
        converged = 1;
        break;
    end
    x = xnew(powsys.bus.busnew);
    iter = iter + 1;
end
Vc = x;
end

