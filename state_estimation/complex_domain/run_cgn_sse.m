function [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, sesettings, varargin)
info.method = 'Gauss-Newton in Complex Variables';

% -------------------- Eq. for the slack bus (no PMU case) ----------------
nopmu = ~meas.num.pmu;
if nopmu
    slackIdx = powsys.num.islack;
else
    slackIdx = [];
end
% -------------------------------------------------------------------------

% --------------------- Setting additional variables ----------------------
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

% ----------------------- Vector of measurement values --------------------
z = [ meas.pmu.m(meas.pmu.IijO) .* exp(1i .* meas.pmu.a(meas.pmu.IijO));
      meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* exp(1i .* meas.pmu.a(meas.pmu.Iinj));
      meas.pmu.m(meas.pmu.IijO) .* exp(-1i .* meas.pmu.a(meas.pmu.IijO));
      meas.pmu.m(meas.pmu.Iij) .* exp(-1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(-1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* exp(-1i .* meas.pmu.a(meas.pmu.Iinj));
      meas.scada.m(meas.scada.pijO);
      meas.scada.m(meas.scada.pij);
      meas.scada.m(meas.scada.qijO);
      meas.scada.m(meas.scada.qij);
      meas.scada.m(meas.scada.pinj);
      meas.scada.m(meas.scada.qinj); 
      meas.scada.m(meas.scada.IijmO);
      meas.scada.m(meas.scada.Iijm);
      meas.scada.m(meas.scada.vm)
      powsys.bus.Vai(slackIdx)
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
if strcmp(sesettings.mweights(1), "pmuscadaratio")
    wIdx = [1:2 * meas.num.pmu + meas.num.scada, (slackIdx/slackIdx) * (2 * meas.num.pmu + meas.num.scada + 1)];
    W = sparse(wIdx, wIdx, [ str2double(sesettings.mweights(2)) .* ones( 2 * meas.num.pmu, 1);...
               ones(meas.num.scada, 1); (slackIdx/slackIdx) * 100 ]);
elseif strcmp(sesettings.mweights(1), "deviceinfo")
    
end
% -------------------------------------------------------------------------

% ------------------------- Construct C11 matrix --------------------------
[ rowInj, colInj ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));  

% ---------------------------- Row indices --------------------------------
iC11IijO = [ 1:meas.num.pIijO, 1:meas.num.pIijO ];
iC11Iij = [ (meas.num.pIijO  + (1:meas.num.pIij)), (meas.num.pIijO  + (1:meas.num.pIij)) ];
iC11V = (meas.num.pIijO + meas.num.pIij + (1:meas.num.pV));
iC11Inj = meas.num.pIijO + meas.num.pIij + meas.num.pV + rowInj;
% -------------------------------------------------------------------------

% ------------------------ Column indices ---------------------------------
jC11Iij0 = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.IijO)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.IijO)); ];
jC11Iij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.Iij));
                 ];
jC11V = meas.pmu.onbus(meas.pmu.v);
jC11Inj = colInj;
% -------------------------------------------------------------------------

% ------------------------- Values of elements ----------------------------
vC11IijO = [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.IijO));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.IijO)); ];
vC11Iij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); 
    ];
vC11V = ones(meas.num.pV, 1);
vC11Inj = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));
% -------------------------------------------------------------------------

C11 = sparse(...
             [ iC11IijO, iC11Iij, iC11V, iC11Inj  ],... 
             [ jC11Iij0; jC11Iij; jC11V; jC11Inj  ], ...
             [ vC11IijO; vC11Iij; vC11V; vC11Inj  ], ...
             meas.num.pmu, powsys.num.bus...
         );
% -------------------------------------------------------------------------

% ------------------------- Construct C12 Matrix -------------------------- 
C12 = sparse(meas.num.pmu, powsys.num.bus);
% -------------------------------------------------------------------------

% -------------------------- Construct D Matrix ---------------------------
[ rowPinj, colPinj ] = find(powsys.ybus.yij(meas.scada.onbus(meas.scada.pinj),:));
[ rowQinj, colQinj ] = find(powsys.ybus.yij(meas.scada.onbus(meas.scada.qinj),:));
accI = 0;
% ---------------------------- Row indices --------------------------------
iDpbrO = [ 1:meas.num.sPijO, 1:meas.num.sPijO ];
accI = accI + meas.num.sPijO;
iDpbr = [ (accI + (1:meas.num.sPij)), (accI + (1:meas.num.sPij)) ];
accI = accI + meas.num.sPij;
iDqbrO = [ (accI + (1:meas.num.sQijO)), (accI + (1:meas.num.sQijO)) ];
accI = accI + meas.num.sQijO;
iDqbr = [ (accI + (1:meas.num.sQij)), (accI + (1:meas.num.sQij)) ];
accI = accI + meas.num.sQij;
iDpinj = [ (accI + (1:meas.num.sPinj)), accI + rowPinj' ];
accI = accI + meas.num.sPinj;
iDqinj = [ (accI + (1:meas.num.sQinj)), accI + rowQinj' ];
accI = accI + meas.num.sQinj;
iDibrMO = [ (accI + (1:meas.num.sIijmO)), (accI + (1:meas.num.sIijmO)) ];
accI = accI + meas.num.sIijmO;
iDibrM = [ (accI + (1:meas.num.sIijm)), (accI + (1:meas.num.sIijm)) ];
accI = accI + meas.num.sIijm;
iDvm = (accI + (1:meas.num.sVm));
if nopmu
    iDsl = meas.num.sVm + accI + 1;
else
    iDsl = [];
end
% -------------------------------------------------------------------------

% --------------------------- Column indices ------------------------------
jDpbrO = [ powsys.branch.i(-meas.scada.loc(meas.scada.pijO));...
           powsys.branch.j(-meas.scada.loc(meas.scada.pijO))];
jDpbr = [ powsys.branch.i(meas.scada.loc(meas.scada.pij));...
          powsys.branch.j(meas.scada.loc(meas.scada.pij))];
jDqbrO = [ powsys.branch.i(-meas.scada.loc(meas.scada.qijO));...
           powsys.branch.j(-meas.scada.loc(meas.scada.qijO))];
jDqbr = [ powsys.branch.i(meas.scada.loc(meas.scada.qij));...
          powsys.branch.j(meas.scada.loc(meas.scada.qij))];   
jDpinj = [ meas.scada.onbus(meas.scada.pinj); colPinj ];
jDqinj = [ meas.scada.onbus(meas.scada.qinj); colQinj ];
jDibrMO = [ powsys.branch.i(-meas.scada.loc(meas.scada.IijmO));...
            powsys.branch.j(-meas.scada.loc(meas.scada.IijmO)) ];
jDibrM = [ powsys.branch.i(meas.scada.loc(meas.scada.Iijm));...
           powsys.branch.j(meas.scada.loc(meas.scada.Iijm)) ];
jDvm = (meas.scada.onbus(meas.scada.vm));

jDsl = slackIdx;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
while iter < sesettings.maxNumberOfIter  
    % ----------------- Branch currents and injected currents -------------
    Ii = powsys.ybus.y * x(powsys.bus.busnew);
    Si = x(powsys.bus.busnew) .* conj(Ii);
    Iij = sum([ powsys.ybus.fromfrom, powsys.ybus.fromto ] .* [ x(powsys.branch.i), x(powsys.branch.j)], 2);
    Sij = x(powsys.branch.i) .* conj(Iij);
    Iji = sum([ powsys.ybus.tofrom, powsys.ybus.toto ] .* [ x(powsys.branch.i), x(powsys.branch.j)], 2);
    Sji = x(powsys.branch.j) .* conj(Iji); 
    % ---------------------------------------------------------------------
    
    % ----------------------- D update & construct  -----------------------
    vDpbrO = [  1/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.pijO)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.pijO)))); ...
                1/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.pijO))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.pijO)))) + ...
                conj(Iji(-meas.scada.loc(meas.scada.pijO))))
                ];
    vDpbr = [  1/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.pij)) ...
               .* x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.pij)))) + ...
               conj(Iij(meas.scada.loc(meas.scada.pij)))); ...
               1/2 .* (powsys.ybus.fromto(meas.scada.loc(meas.scada.pij))) .* ...
               x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.pij))))
                ];  
    vDqbrO = [  1i/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.qijO)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.qijO)))); ...
                1i/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.qijO))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.qijO)))) - ...
                conj(Iji(-meas.scada.loc(meas.scada.qijO))))
                ];
    vDqbr = [  1i/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.qij)) ...
               .* x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.qij)))) - ...
               conj(Iij(meas.scada.loc(meas.scada.qij)))); ...
               1i/2 .* (powsys.ybus.fromto(meas.scada.loc(meas.scada.qij))) .* ...
               x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.qij))))
                ]; 
       
    vDpinj = [ 1/2 .* (conj(Ii(meas.scada.onbus(meas.scada.pinj))) + ...
               powsys.ybus.ydiag(meas.scada.onbus(meas.scada.pinj)) .* ...
               x(powsys.num.bus + meas.scada.onbus(meas.scada.pinj)));...
                nonzeros(1/2 .*  conj(x(meas.scada.onbus(meas.scada.pinj))) .* ...
                powsys.ybus.yij(meas.scada.onbus(meas.scada.pinj), :))];
    vDqinj = [ 1i/2 .* (-conj(Ii(meas.scada.onbus(meas.scada.qinj))) + ...
               powsys.ybus.ydiag(meas.scada.onbus(meas.scada.qinj)) .* ...
               x(powsys.num.bus + meas.scada.onbus(meas.scada.qinj)));...
               nonzeros(1i/2 .*  conj(x(meas.scada.onbus(meas.scada.qinj))) .* ...
                powsys.ybus.yij(meas.scada.onbus(meas.scada.qinj), :))];     
    vDibrMO = [ 1/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.IijmO)) ...
                .* exp(-1i .* angle(Iji(-meas.scada.loc(meas.scada.IijmO))));
                1/2 .* powsys.ybus.toto(-meas.scada.loc(meas.scada.IijmO)) ...
                .* exp(-1i .* angle(Iji(-meas.scada.loc(meas.scada.IijmO)))); ];
    vDibrM =  [ 1/2 .* powsys.ybus.fromfrom(meas.scada.loc(meas.scada.Iijm)) ...
                .* exp(-1i .* angle(Iij(meas.scada.loc(meas.scada.Iijm))));
                1/2 .* powsys.ybus.fromto(meas.scada.loc(meas.scada.Iijm)) ...
                .* exp(-1i .* angle(Iij(meas.scada.loc(meas.scada.Iijm)))); ];
    vDvm = 1/2 .* exp(-1i .* angle(meas.scada.onbus(meas.scada.vm)));
    if nopmu
        vDsl = -1i/2;
    else
        vDsl = [];
    end
    
    D = sparse(...
               [ iDpbrO, iDpbr, iDqbrO, iDqbr, iDpinj, iDqinj, iDibrMO, iDibrM, iDvm, iDsl ], ...
               [ jDpbrO; jDpbr; jDqbrO; jDqbr; jDpinj; jDqinj; jDibrMO; jDibrM; jDvm; jDsl ], ...
               [ vDpbrO; vDpbr; vDqbrO; vDqbr; vDpinj; vDqinj; vDibrMO; vDibrM; vDvm; vDsl ] );
    % ---------------------------------------------------------------------
    
    % --------------------------- h <-> PMU -------------------------------
    c = [   Iji(-meas.pmu.loc(meas.pmu.IijO));
            Iij(meas.pmu.loc(meas.pmu.Iij));
            x(meas.pmu.onbus(meas.pmu.v));
            Ii(meas.pmu.onbus(meas.pmu.Iinj))
          ];
    % ---------------------------------------------------------------------
    
    % ----------------------------- d <-> SCADA ---------------------------
    d = [ real(Sji(-meas.scada.loc(meas.scada.pijO)));
          real(Sij(meas.scada.loc(meas.scada.pij)));
          imag(Sji(-meas.scada.loc(meas.scada.qijO)));
          imag(Sij(meas.scada.loc(meas.scada.qij)));
          real(Si(meas.scada.onbus(meas.scada.pinj)));
          imag(Si(meas.scada.onbus(meas.scada.qinj)));
          abs(Iji(-meas.scada.loc(meas.scada.IijmO)));
          abs(Iij(meas.scada.loc(meas.scada.Iijm)));
          abs(x(meas.scada.onbus(meas.scada.vm)))
          angle(x(slackIdx))
          ];
    % ---------------------------------------------------------------------
    H = [  C11         C12
         conj(C12)   conj(C11)
           D      conj(D)];
    h = [ c; conj(c); d ];
    
    % Solve
    r = z - h;
    dx = (H' * W * H) \ (H' * W * r);
    x = x + dx; 
    if max(abs(dx)) < sesettings.eps
        converged = 1;
        break
    else
        iter = iter + 1;
    end
end
% ----------------------- Estimator results -------------------------------
Vc = x;
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end

