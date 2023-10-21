function [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, sesettings, varargin)
info.method = 'Gauss-Newton in Complex Variables';

% -------------------- Eq. for the slack bus (no PMU case) ----------------
nopmu = ~meas.num.pmu;
noscada = ~meas.num.scada;
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
z = [ meas.pmu.m(meas.pmu.Iji) .* exp(1i .* meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* exp(1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* exp(1i .* meas.pmu.a(meas.pmu.Iinj));
      meas.pmu.m(meas.pmu.Iji) .* exp(-1i .* meas.pmu.a(meas.pmu.Iji));
      meas.pmu.m(meas.pmu.Iij) .* exp(-1i .* meas.pmu.a(meas.pmu.Iij));
      meas.pmu.m(meas.pmu.v) .* exp(-1i .* meas.pmu.a(meas.pmu.v));
      meas.pmu.m(meas.pmu.Iinj) .* exp(-1i .* meas.pmu.a(meas.pmu.Iinj));
      zeros(2 * powsys.num.zi * sesettings.virtual, 1);
      meas.scada.m(meas.scada.pji);
      meas.scada.m(meas.scada.pij);
      meas.scada.m(meas.scada.qji);
      meas.scada.m(meas.scada.qij);
      meas.scada.m(meas.scada.pinj);
      meas.scada.m(meas.scada.qinj); 
      meas.scada.m(meas.scada.Ijim);
      meas.scada.m(meas.scada.Iijm);
      meas.scada.m(meas.scada.vm)
      powsys.bus.Vai(slackIdx)
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
if strcmp(sesettings.mweights(1), "pmuscadaratio")
    wIdx = [1:2 * meas.num.pmu + meas.num.scada + sesettings.virtual * 2 * powsys.num.zi,...
           (slackIdx/slackIdx) * (2 * meas.num.pmu + meas.num.scada + 2 * sesettings.virtual * powsys.num.zi + 1)];
    W = sparse(wIdx, wIdx, [ str2double(sesettings.mweights(2)) .* ones( 2 * meas.num.pmu, 1);...
               25 .* ones(2 * powsys.num.zi * sesettings.virtual, 1); ones(meas.num.scada, 1); (slackIdx/slackIdx) * 25 ]);
elseif strcmp(sesettings.mweights(1), "deviceinfo")
    W = [];
end
% -------------------------------------------------------------------------

% ------------------------- Construct C11 matrix --------------------------
[ rowInj, colInj ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));  
% ---------------------------- Row indices --------------------------------
iC11Iji = [ 1:meas.num.pIji, 1:meas.num.pIji ];
iC11Iij = [ (meas.num.pIji  + (1:meas.num.pIij)), (meas.num.pIji  + (1:meas.num.pIij)) ];
iC11V = (meas.num.pIji + meas.num.pIij + (1:meas.num.pV));
iC11Inj = meas.num.pIji + meas.num.pIij + meas.num.pV + rowInj;
% -------------------------------------------------------------------------

% ------------------------ Column indices ---------------------------------
jC11Iji = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.Iji)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.Iji)); ];
jC11Iij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.Iij)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.Iij));
                 ];
jC11V = meas.pmu.onbus(meas.pmu.v);
jC11Inj = colInj;
% -------------------------------------------------------------------------

% ------------------------- Values of elements ----------------------------
vC11Iji = [  powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji));...
             powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)); 
             ];
vC11Iij =  [ powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
             powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); ];
vC11V = ones(meas.num.pV, 1);
vC11Inj = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));
% -------------------------------------------------------------------------

C11 = sparse(...
             [ iC11Iji, iC11Iij, iC11V, iC11Inj  ],... 
             [ jC11Iji; jC11Iij; jC11V; jC11Inj  ], ...
             [ vC11Iji; vC11Iij; vC11V; vC11Inj  ], ...
             meas.num.pmu, powsys.num.bus );
% -------------------------------------------------------------------------

% ---------------- OPTIONAL - Virtual current injection -------------------
if sesettings.virtual
    J = powsys.ybus.y(powsys.bus.zi,:); 
    Zj = sparse(powsys.num.zi, powsys.num.bus);
else
    J = [];
    Zj = [];
end
% -------------------------------------------------------------------------
% ------------------------- Construct C12 Matrix -------------------------- 
C12 = sparse(meas.num.pmu, powsys.num.bus);
% -------------------------------------------------------------------------

% -------------------------- Construct D Matrix ---------------------------
[ rowPinj, colPinj ] = find(powsys.ybus.yij(meas.scada.onbus(meas.scada.pinj),:));
[ rowQinj, colQinj ] = find(powsys.ybus.yij(meas.scada.onbus(meas.scada.qinj),:));
accI = 0;
% ---------------------------- Row indices --------------------------------
iDpji = [ 1:meas.num.sPji, 1:meas.num.sPji ];
accI = accI + meas.num.sPji;
iDpij = [ (accI + (1:meas.num.sPij)), (accI + (1:meas.num.sPij)) ];
accI = accI + meas.num.sPij;
iDqji = [ (accI + (1:meas.num.sQji)), (accI + (1:meas.num.sQji)) ];
accI = accI + meas.num.sQji;
iDqij = [ (accI + (1:meas.num.sQij)), (accI + (1:meas.num.sQij)) ];
accI = accI + meas.num.sQij;
iDpinj = [ (accI + (1:meas.num.sPinj)), accI + rowPinj' ];
accI = accI + meas.num.sPinj;
iDqinj = [ (accI + (1:meas.num.sQinj)), accI + rowQinj' ];
accI = accI + meas.num.sQinj;
iDIjim = [ (accI + (1:meas.num.sIjim)), (accI + (1:meas.num.sIjim)) ];
accI = accI + meas.num.sIjim;
iDIijm = [ (accI + (1:meas.num.sIijm)), (accI + (1:meas.num.sIijm)) ];
accI = accI + meas.num.sIijm;
iDvm = (accI + (1:meas.num.sVm));
if nopmu
    iDsl = meas.num.sVm + accI + 1;
else
    iDsl = [];
end
% -------------------------------------------------------------------------

% --------------------------- Column indices ------------------------------
jDpji = [ powsys.branch.i(-meas.scada.loc(meas.scada.pji));...
           powsys.branch.j(-meas.scada.loc(meas.scada.pji))];
jDpij = [ powsys.branch.i(meas.scada.loc(meas.scada.pij));...
          powsys.branch.j(meas.scada.loc(meas.scada.pij))];
jDqji = [ powsys.branch.i(-meas.scada.loc(meas.scada.qji));...
           powsys.branch.j(-meas.scada.loc(meas.scada.qji))];
jDqij = [ powsys.branch.i(meas.scada.loc(meas.scada.qij));...
          powsys.branch.j(meas.scada.loc(meas.scada.qij))];   
jDpinj = [ meas.scada.onbus(meas.scada.pinj); colPinj ];
jDqinj = [ meas.scada.onbus(meas.scada.qinj); colQinj ];
jDIjim = [ powsys.branch.i(-meas.scada.loc(meas.scada.Ijim));...
            powsys.branch.j(-meas.scada.loc(meas.scada.Ijim)) ];
jDIijm = [ powsys.branch.i(meas.scada.loc(meas.scada.Iijm));...
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
    vDpji = [  1/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.pji)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.pji)))); ...
                1/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.pji))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.pji)))) + ...
                conj(Iji(-meas.scada.loc(meas.scada.pji))))
                ];
    vDpij = [  1/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.pij)) ...
               .* x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.pij)))) + ...
               conj(Iij(meas.scada.loc(meas.scada.pij)))); ...
               1/2 .* (powsys.ybus.fromto(meas.scada.loc(meas.scada.pij))) .* ...
               x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.pij))))
                ];  
    vDqji = [  1i/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.qji)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.qji)))); ...
                1i/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.qji))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.qji)))) - ...
                conj(Iji(-meas.scada.loc(meas.scada.qji))))
                ];
    vDqij = [  1i/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.qij)) ...
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
    vDIjim = [ 1/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.Ijim)) ...
                .* exp(-1i .* angle(Iji(-meas.scada.loc(meas.scada.Ijim))));
                1/2 .* powsys.ybus.toto(-meas.scada.loc(meas.scada.Ijim)) ...
                .* exp(-1i .* angle(Iji(-meas.scada.loc(meas.scada.Ijim)))); ];
    vDIijm =  [ 1/2 .* powsys.ybus.fromfrom(meas.scada.loc(meas.scada.Iijm)) ...
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
               [ iDpij, iDpji, iDqji, iDqij, iDpinj, iDqinj, iDIijm, iDIjim, iDvm, iDsl ], ...
               [ jDpij, jDpji, jDqji, jDqij, jDpinj, jDqinj, jDIijm, jDIjim, jDvm, jDsl ], ...
               [ vDpij, vDpji, vDqji, vDqij, vDpinj, vDqinj, vDIijm, vDIjim, vDvm, vDsl ] );
    % ---------------------------------------------------------------------
    
    % --------------------------- h <-> PMU -------------------------------
    c = [   Iji(-meas.pmu.loc(meas.pmu.Iji));
            Iij(meas.pmu.loc(meas.pmu.Iij));
            x(meas.pmu.onbus(meas.pmu.v));
            Ii(meas.pmu.onbus(meas.pmu.Iinj))
          ];
    % ---------------------------------------------------------------------
    
    % ---------------------------- j <-> ZI -------------------------------
    if sesettings.virtual
        j = Ii(powsys.bus.zi);
    else
        j = [];
    end
    % ---------------------------------------------------------------------
    
    % ----------------------------- d <-> SCADA ---------------------------
    d = [ real(Sji(-meas.scada.loc(meas.scada.pji)));
          real(Sij(meas.scada.loc(meas.scada.pij)));
          imag(Sji(-meas.scada.loc(meas.scada.qji)));
          imag(Sij(meas.scada.loc(meas.scada.qij)));
          real(Si(meas.scada.onbus(meas.scada.pinj)));
          imag(Si(meas.scada.onbus(meas.scada.qinj)));
          abs(Iji(-meas.scada.loc(meas.scada.Ijim)));
          abs(Iij(meas.scada.loc(meas.scada.Iijm)));
          abs(x(meas.scada.onbus(meas.scada.vm)))
          angle(x(slackIdx))
          ];
    % ---------------------------------------------------------------------
    H = [  C11         C12
         conj(C12)   conj(C11)
            J          Zj
            Zj       conj(J)
           D      conj(D)];
    h = [ c; conj(c); j; conj(j); d ];
    
    % Solve
    r = z - h;
    dx = (H' * W * H) \ (H' * W * r);
    x = x + dx; 
    if max(abs(dx)) < sesettings.eps || noscada
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

