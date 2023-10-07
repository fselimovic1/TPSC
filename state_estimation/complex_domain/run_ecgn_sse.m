function [ Vc, iter, converged, info ] = run_ecgn_sse(powsys, meas, sesettings, varargin)
info.method = 'Gauss-Newton (full version) in Complex Variables';

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
      meas.scada.m(meas.scada.sji.p) + 1i * meas.scada.m(meas.scada.sji.q)
      meas.scada.m(meas.scada.sij.p) + 1i * meas.scada.m(meas.scada.sij.q)
      meas.scada.m(meas.scada.sinj.p) + 1i * meas.scada.m(meas.scada.sinj.q)
      meas.scada.m(meas.scada.sji.p) - 1i * meas.scada.m(meas.scada.sji.q)
      meas.scada.m(meas.scada.sij.p) - 1i * meas.scada.m(meas.scada.sij.q)
      meas.scada.m(meas.scada.sinj.p) - 1i * meas.scada.m(meas.scada.sinj.q)
      meas.scada.m(meas.scada.oddPji);
      meas.scada.m(meas.scada.oddPij);
      meas.scada.m(meas.scada.oddQji);
      meas.scada.m(meas.scada.oddQij);
      meas.scada.m(meas.scada.oddPinj);
      meas.scada.m(meas.scada.oddQinj); 
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
               2 * str2double(sesettings.mweights(2))  .* ones(2 * powsys.num.zi * sesettings.virtual, 1);...
               ones(meas.num.scada, 1); (slackIdx/slackIdx) * 2 * str2double(sesettings.mweights(2))]);
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
vC11Iji = [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.Iji));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.Iji)); ];
vC11Iij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.Iij));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.Iij)); 
    ];
vC11V = ones(meas.num.pV, 1);
vC11Inj = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.Iinj), :));
% -------------------------------------------------------------------------

C11 = sparse(...
             [ iC11Iji, iC11Iij, iC11V, iC11Inj  ],... 
             [ jC11Iji; jC11Iij; jC11V; jC11Inj  ], ...
             [ vC11Iji; vC11Iij; vC11V; vC11Inj  ], meas.num.pmu, powsys.num.bus );
% -------------------------------------------------------------------------

% ------------------------- Construct C12 Matrix -------------------------- 
C12 = sparse(meas.num.pmu, powsys.num.bus);
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

% --------------------- Complex SCADA measurements (Cs) -------------------
[ rowInj, colInj ] = find(powsys.ybus.y(meas.scada.loc(meas.scada.sinj.p), :));  
% ------------------------- Row indices (Cs11) ----------------------------
iCs11Sji = 1:meas.num.sSji;
iCs11Sij = meas.num.sSji + (1:meas.num.sSij);
iCs11Sinj = meas.num.sSij + meas.num.sSji + (1:meas.num.sSinj);
% -------------------------------------------------------------------------

% ------------------------- Column indices (Cs11) -------------------------
jCs11Sji = powsys.branch.i(-meas.scada.loc(meas.scada.sji.p));
jCs11Sij = powsys.branch.i(meas.scada.loc(meas.scada.sij.p));
jCs11Sinj = meas.scada.loc(meas.scada.sinj.p);
% -------------------------------------------------------------------------

% ------------------------- Row indices (Cs12) ----------------------------
iCs12Sji = [ 1:meas.num.sSji, 1:meas.num.sSji ];
iCs12Sij = meas.num.sSji + [ 1:meas.num.sSij, 1:meas.num.sSij ];
iCs12Sinj = meas.num.sSij + meas.num.sSji + rowInj;
% -------------------------------------------------------------------------

% ------------------------- Column indices (Cs12) -------------------------
jCs12Sji = [ powsys.branch.i(-meas.scada.loc(meas.scada.sji.p));
             powsys.branch.j(-meas.scada.loc(meas.scada.sji.p));
             ];
jCs12Sij = [ powsys.branch.i(meas.scada.loc(meas.scada.sij.p));
             powsys.branch.j(meas.scada.loc(meas.scada.sij.p));
             ];
jCs12Sinj = colInj;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------


% -------------------------- Construct D Matrix ---------------------------
[ rowPinj, colPinj ] = find(powsys.ybus.yij(meas.scada.onbus(meas.scada.oddPinj),:));
[ rowQinj, colQinj ] = find(powsys.ybus.yij(meas.scada.onbus(meas.scada.oddQinj),:));
accI = 0;
% ---------------------------- Row indices --------------------------------
iDpji = [ 1:meas.num.soPji, 1:meas.num.soPji ];
accI = accI + meas.num.soPji;
iDpij = [ (accI + (1:meas.num.soPij)), (accI + (1:meas.num.soPij)) ];
accI = accI + meas.num.soPij;
iDqji = [ (accI + (1:meas.num.soQji)), (accI + (1:meas.num.soQji)) ];
accI = accI + meas.num.soQji;
iDqij = [ (accI + (1:meas.num.soQij)), (accI + (1:meas.num.soQij)) ];
accI = accI + meas.num.soQij;
iDpinj = [ (accI + (1:meas.num.soPinj)), accI + rowPinj' ];
accI = accI + meas.num.soPinj;
iDqinj = [ (accI + (1:meas.num.soQinj)), accI + rowQinj' ];
accI = accI + meas.num.soQinj;
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
jDpji = [ powsys.branch.i(-meas.scada.loc(meas.scada.oddPji));...
           powsys.branch.j(-meas.scada.loc(meas.scada.oddPji))];
jDpij = [ powsys.branch.i(meas.scada.loc(meas.scada.oddPij));...
          powsys.branch.j(meas.scada.loc(meas.scada.oddPij))];
jDqji = [ powsys.branch.i(-meas.scada.loc(meas.scada.oddQji));...
           powsys.branch.j(-meas.scada.loc(meas.scada.oddQji))];
jDqij = [ powsys.branch.i(meas.scada.loc(meas.scada.oddQij));...
          powsys.branch.j(meas.scada.loc(meas.scada.oddQij))];   
jDpinj = [ meas.scada.onbus(meas.scada.oddPinj); colPinj ];
jDqinj = [ meas.scada.onbus(meas.scada.oddQinj); colQinj ];
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
    % ------------------------ Cs11 & Cs12 update -------------------------
    % ---------------------------- Values (Cs11) --------------------------
    vCs11Sji = conj(Iji(-meas.scada.loc(meas.scada.sji.p))); 
    vCs11Sij = conj(Iij(meas.scada.loc(meas.scada.sij.p))); 
    vCs11Sinj = conj(Ii(meas.scada.loc(meas.scada.sinj.p)));
    % ---------------------------------------------------------------------
    Cs11 = sparse([iCs11Sji, iCs11Sij, iCs11Sinj], [jCs11Sji, jCs11Sij, jCs11Sinj] ,...
        [vCs11Sji, vCs11Sij, vCs11Sinj], meas.num.sSji + meas.num.sSij + ...
         meas.num.sSinj + meas.num.sIijm + meas.num.sIjim + meas.num.sVm, powsys.num.bus);

    % ---------------------------- Values (Cs12) --------------------------
    vCs12Sji = [ x(powsys.branch.j(-meas.scada.loc(meas.scada.sji.p))) .* ...
                 conj(powsys.ybus.tofrom(-meas.scada.loc(meas.scada.sji.p)));
                 x(powsys.branch.j(-meas.scada.loc(meas.scada.sji.p))) .* ...
                 conj(powsys.ybus.toto(-meas.scada.loc(meas.scada.sji.p)));
                 ];
    vCs12Sij = [ x(powsys.branch.i(meas.scada.loc(meas.scada.sij.p))) .* ...
                 conj(powsys.ybus.fromfrom(meas.scada.loc(meas.scada.sij.p)));
                 x(powsys.branch.i(meas.scada.loc(meas.scada.sij.p))) .* ...
                 conj(powsys.ybus.fromto(meas.scada.loc(meas.scada.sij.p)));
                 ];
    vCs12Sinj = nonzeros(powsys.ybus.y(meas.scada.loc(meas.scada.sinj.p), :)) .* ...
                x(meas.scada.loc(meas.scada.sinj.p));
    % -------------------------------------------------------------------------
     Cs12 = sparse([iCs12Sji, iCs12Sij, iCs12Sinj], [jCs12Sji, jCs12Sij, jCs12Sinj] ,...
        [vCs12Sji, vCs12Sij, vCs12Sinj], meas.num.sSji + meas.num.sSij + ...
         meas.num.sSinj + meas.num.sIijm + meas.num.sIjim + meas.num.sVm, powsys.num.bus);
    % ---------------------------------------------------------------------
    
    % ----------------------- D update & construct  -----------------------
    vDpji = [  1/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.oddPji)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.oddPji)))); ...
                1/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.oddPji))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.oddPji)))) + ...
                conj(Iji(-meas.scada.loc(meas.scada.oddPji))))
                ];
    vDpij = [  1/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.oddPij)) ...
               .* x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.oddPij)))) + ...
               conj(Iij(meas.scada.loc(meas.scada.oddPij)))); ...
               1/2 .* (powsys.ybus.fromto(meas.scada.loc(meas.scada.oddPij))) .* ...
               x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.oddPij))))
                ];  
    vDqji = [  1i/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.oddQji)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.oddQji)))); ...
                1i/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.oddQji))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.oddQji)))) - ...
                conj(Iji(-meas.scada.loc(meas.scada.oddQji))))
                ];
    vDqij = [  1i/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.oddQij)) ...
               .* x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.oddQij)))) - ...
               conj(Iij(meas.scada.loc(meas.scada.oddQij)))); ...
               1i/2 .* (powsys.ybus.fromto(meas.scada.loc(meas.scada.oddQij))) .* ...
               x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.oddQij))))
                ]; 
       
    vDpinj = [ 1/2 .* (conj(Ii(meas.scada.onbus(meas.scada.oddPinj))) + ...
               powsys.ybus.ydiag(meas.scada.onbus(meas.scada.oddPinj)) .* ...
               x(powsys.num.bus + meas.scada.onbus(meas.scada.oddPinj)));...
                nonzeros(1/2 .*  conj(x(meas.scada.onbus(meas.scada.oddPinj))) .* ...
                powsys.ybus.yij(meas.scada.onbus(meas.scada.oddPinj), :))];
    vDqinj = [ 1i/2 .* (-conj(Ii(meas.scada.onbus(meas.scada.oddQinj))) + ...
               powsys.ybus.ydiag(meas.scada.onbus(meas.scada.oddQinj)) .* ...
               x(powsys.num.bus + meas.scada.onbus(meas.scada.oddQinj)));...
               nonzeros(1i/2 .*  conj(x(meas.scada.onbus(meas.scada.oddQinj))) .* ...
                powsys.ybus.yij(meas.scada.onbus(meas.scada.oddQinj), :))];     
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
    cs = [ Sji(-meas.scada.loc(meas.scada.sji.p))
           Sij(meas.scada.loc(meas.scada.sij.q)) 
           Si(meas.scada.loc(meas.scada.sinj.q))
    ];
    d = [ real(Sji(-meas.scada.loc(meas.scada.oddPji)));
          real(Sij(meas.scada.loc(meas.scada.oddPij)));
          imag(Sji(-meas.scada.loc(meas.scada.oddQji)));
          imag(Sij(meas.scada.loc(meas.scada.oddQij)));
          real(Si(meas.scada.onbus(meas.scada.oddPinj)));
          imag(Si(meas.scada.onbus(meas.scada.oddQinj)));
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
            Cs11         Cs12
         conj(Cs12)   conj(Cs11)
           D      conj(D)];
    h = [ c; conj(c); j; conj(j); cs; conj(cs); d ];
    
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

