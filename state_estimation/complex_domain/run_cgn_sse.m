function [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, sesettings, varargin)
info.method = 'Gauss-Newton in Complex Variables';
% --------------------- Setting additional variables ----------------------
iter = 1;
converged = 0;
% -------------------------------------------------------------------------

% --------------------- Initialize state variables ------------------------
if nargin == 4
    x0 = varargin{1};
    x = [ x0;
        conj(x0) 
        ];
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
z = [ meas.pmu.m(meas.pmu.ibranchO) .* exp(1i .* meas.pmu.a(meas.pmu.ibranchO));
      meas.pmu.m(meas.pmu.ibranch) .* exp(1i .* meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* exp(1i .* meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.inj) .* exp(1i .* meas.pmu.a(meas.pmu.inj));
      meas.pmu.m(meas.pmu.ibranchO) .* exp(-1i .* meas.pmu.a(meas.pmu.ibranchO));
      meas.pmu.m(meas.pmu.ibranch) .* exp(-1i .* meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* exp(-1i .* meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.inj) .* exp(-1i .* meas.pmu.a(meas.pmu.inj));
      meas.scada.m(meas.scada.pbranchO);
      meas.scada.m(meas.scada.pbranch);
      meas.scada.m(meas.scada.qbranchO);
      meas.scada.m(meas.scada.qbranch);
      meas.scada.m(meas.scada.pinj);
      meas.scada.m(meas.scada.qinj); 
      meas.scada.m(meas.scada.ibrMO);
      meas.scada.m(meas.scada.ibrM);
      meas.scada.m(meas.scada.vm)
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
W = sparse(1:2 * meas.num.pmu + meas.num.scada, 1:2 * meas.num.pmu + meas.num.scada, ...
           [ str2double(sesettings.mweights(2)) .* ones( 2 * meas.num.pmu, 1);...
           ones(meas.num.scada, 1) ]);
% -------------------------------------------------------------------------

% ------------------------- Construct C11 matrix --------------------------
[ rowInj, colInj ] = find(powsys.ybus.y(meas.pmu.loc(meas.pmu.inj), :));  

% ---------------------------- Row indices --------------------------------
iC11IijO = [ 1:meas.num.pIijO, 1:meas.num.pIijO ];
iC11Iij = [ (meas.num.pIijO  + (1:meas.num.pIij)), (meas.num.pIijO  + (1:meas.num.pIij)) ];
iC11V = (meas.num.pIijO + meas.num.pIij + (1:meas.num.pV));
iC11Inj = meas.num.pIijO + meas.num.pIij + meas.num.pV + rowInj;
% -------------------------------------------------------------------------

% ------------------------ Column indices ---------------------------------
jC11Iij0 = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.ibranchO)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.ibranchO)); ];
jC11Iij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.ibranch)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.ibranch));
                 ];
jC11V = meas.pmu.onbus(meas.pmu.vnode);
jC11Inj = colInj;
% -------------------------------------------------------------------------

% ------------------------- Values of elements ----------------------------
vC11IijO = [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchO));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchO)); ];
vC11Iij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)); 
    ];
vC11V = ones(meas.num.pV, 1);
vC11Inj = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.inj), :));
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
iDpbrO = [ 1:meas.num.sPbr0, 1:meas.num.sPbr0 ];
accI = accI + meas.num.sPbr0;
iDpbr = [ (accI + (1:meas.num.sPbr)), (accI + (1:meas.num.sPbr)) ];
accI = accI + meas.num.sPbr;
iDqbrO = [ (accI + (1:meas.num.sQbr0)), (accI + (1:meas.num.sQbr0)) ];
accI = accI + meas.num.sQbr0;
iDqbr = [ (accI + (1:meas.num.sQbr)), (accI + (1:meas.num.sQbr)) ];
accI = accI + meas.num.sQbr;
iDpinj = [ (accI + (1:meas.num.sPinj)), accI + rowPinj' ];
accI = accI + meas.num.sPinj;
iDqinj = [ (accI + (1:meas.num.sQinj)), accI + rowQinj' ];
accI = accI + meas.num.sQinj;
iDibrMO = [ (accI + (1:meas.num.sIbrM0)), (accI + (1:meas.num.sIbrM0)) ];
accI = accI + meas.num.sIbrM0;
iDibrM = [ (accI + (1:meas.num.sIbrM)), (accI + (1:meas.num.sIbrM)) ];
accI = accI + meas.num.sIbrM;
iDvm = (accI + (1:meas.num.sVm));
% -------------------------------------------------------------------------

% --------------------------- Column indices ------------------------------
jDpbrO = [ powsys.branch.i(-meas.scada.loc(meas.scada.pbranchO));...
           powsys.branch.j(-meas.scada.loc(meas.scada.pbranchO))];
jDpbr = [ powsys.branch.i(meas.scada.loc(meas.scada.pbranch));...
          powsys.branch.j(meas.scada.loc(meas.scada.pbranch))];
jDqbrO = [ powsys.branch.i(-meas.scada.loc(meas.scada.qbranchO));...
           powsys.branch.j(-meas.scada.loc(meas.scada.qbranchO))];
jDqbr = [ powsys.branch.i(meas.scada.loc(meas.scada.qbranch));...
          powsys.branch.j(meas.scada.loc(meas.scada.qbranch))];   
jDpinj = [ meas.scada.onbus(meas.scada.pinj); colPinj ];
jDqinj = [ meas.scada.onbus(meas.scada.qinj); colQinj ];
jDibrMO = [ powsys.branch.i(-meas.scada.loc(meas.scada.ibrMO));...
            powsys.branch.j(-meas.scada.loc(meas.scada.ibrMO)) ];
jDibrM = [ powsys.branch.i(meas.scada.loc(meas.scada.ibrM));...
           powsys.branch.j(meas.scada.loc(meas.scada.ibrM)) ];
jDvm = (meas.scada.onbus(meas.scada.vm));
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
    vDpbrO = [  1/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.pbranchO)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.pbranchO)))); ...
                1/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.pbranchO))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.pbranchO)))) + ...
                conj(Iji(-meas.scada.loc(meas.scada.pbranchO))))
                ];
    vDpbr = [  1/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.pbranch)) ...
               .* x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.pbranch)))) + ...
               conj(Iij(meas.scada.loc(meas.scada.pbranch)))); ...
               1/2 .* (powsys.ybus.fromto(meas.scada.loc(meas.scada.pbranch))) .* ...
               x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.pbranch))))
                ];  
    vDqbrO = [  1i/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.qbranchO)) ...
               .* x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.qbranchO)))); ...
                1i/2 .* ((powsys.ybus.toto(-meas.scada.loc(meas.scada.qbranchO))) .* ...
                x(powsys.num.bus + (powsys.branch.j(-meas.scada.loc(meas.scada.qbranchO)))) - ...
                conj(Iji(-meas.scada.loc(meas.scada.qbranchO))))
                ];
    vDqbr = [  1i/2 .* (powsys.ybus.fromfrom(meas.scada.loc(meas.scada.qbranch)) ...
               .* x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.qbranch)))) - ...
               conj(Iij(meas.scada.loc(meas.scada.qbranch)))); ...
               1i/2 .* (powsys.ybus.fromto(meas.scada.loc(meas.scada.qbranch))) .* ...
               x(powsys.num.bus + (powsys.branch.i(meas.scada.loc(meas.scada.qbranch))))
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
    vDibrMO = [ 1/2 .* powsys.ybus.tofrom(-meas.scada.loc(meas.scada.ibrMO)) ...
                .* exp(-1i .* angle(Iji(-meas.scada.loc(meas.scada.ibrMO))));
                1/2 .* powsys.ybus.toto(-meas.scada.loc(meas.scada.ibrMO)) ...
                .* exp(-1i .* angle(Iji(-meas.scada.loc(meas.scada.ibrMO)))); ];
    vDibrM =  [ 1/2 .* powsys.ybus.fromfrom(meas.scada.loc(meas.scada.ibrM)) ...
                .* exp(-1i .* angle(Iij(meas.scada.loc(meas.scada.ibrM))));
                1/2 .* powsys.ybus.fromto(meas.scada.loc(meas.scada.ibrM)) ...
                .* exp(-1i .* angle(Iij(meas.scada.loc(meas.scada.ibrM)))); ];
    vDvm = 1/2 .* exp(-1i .* angle(meas.scada.onbus(meas.scada.vm)));
    
    D = sparse(...
               [ iDpbrO, iDpbr, iDqbrO, iDqbr, iDpinj, iDqinj, iDibrMO, iDibrM, iDvm ], ...
               [ jDpbrO; jDpbr; jDqbrO; jDqbr; jDpinj; jDqinj; jDibrMO; jDibrM; jDvm ], ...
               [ vDpbrO; vDpbr; vDqbrO; vDqbr; vDpinj; vDqinj; vDibrMO; vDibrM; vDvm ] );
    % ---------------------------------------------------------------------
    
    % --------------------------- h <-> PMU -------------------------------
    c = [   Iji(-meas.pmu.loc(meas.pmu.ibranchO));
            Iij(meas.pmu.loc(meas.pmu.ibranch));
            x(meas.pmu.onbus(meas.pmu.vnode));
            Ii(meas.pmu.onbus(meas.pmu.inj))
          ];
    % ---------------------------------------------------------------------
    
    % ----------------------------- d <-> SCADA ---------------------------
    d = [ real(Sji(-meas.scada.loc(meas.scada.pbranchO)));
          real(Sij(meas.scada.loc(meas.scada.pbranch)));
          imag(Sji(-meas.scada.loc(meas.scada.qbranchO)));
          imag(Sij(meas.scada.loc(meas.scada.qbranch)));
          real(Si(meas.scada.onbus(meas.scada.pinj)));
          imag(Si(meas.scada.onbus(meas.scada.qinj)));
          abs(Iji(-meas.scada.loc(meas.scada.ibrMO)));
          abs(Iij(meas.scada.loc(meas.scada.ibrM)));
          abs(x(meas.scada.onbus(meas.scada.vm)))
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
Vc = x(powsys.bus.busnew);
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end

