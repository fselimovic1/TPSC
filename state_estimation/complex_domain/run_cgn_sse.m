function [ Vc, iter, converged, info ] = run_cgn_sse(powsys, meas, sesettings)
info.method = 'Gauss-Newton in Complex Variables';
% --------------------- Setting additional variables ----------------------
iter = 1;
converged = 0;
% -------------------------------------------------------------------------

% --------------------- Initialize state variables ------------------------
if sesettings.flatStart
    x = ones( 2 * powsys.num.bus, 1);
else
    x = [ powsys.bus.Vmi .* (1i * powsys.bus.Vai);
          powsys.bus.Vmi .* (-1i * powsys.bus.Vai);
         ];
end
% -------------------------------------------------------------------------

% ----------------------- Vector of measurement values --------------------
z = [ meas.pmu.m(meas.pmu.ibranchOpp) .* exp(1i .* meas.pmu.a(meas.pmu.ibranchOpp));
      meas.pmu.m(meas.pmu.ibranch) .* exp(1i .* meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* exp(1i .* meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.inj) .* exp(1i .* meas.pmu.a(meas.pmu.inj));
      meas.pmu.m(meas.pmu.ibranchOpp) .* exp(-1i .* meas.pmu.a(meas.pmu.ibranchOpp));
      meas.pmu.m(meas.pmu.ibranch) .* exp(-1i .* meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* exp(-1i .* meas.pmu.a(meas.pmu.vnode));
      meas.pmu.m(meas.pmu.inj) .* exp(-1i .* meas.pmu.a(meas.pmu.inj));
      meas.scada.m
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
jC11Iij0 = [     powsys.branch.i(-meas.pmu.loc(meas.pmu.ibranchOpp)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.ibranchOpp)); ];
jC11Iij =  [     powsys.branch.i(meas.pmu.loc(meas.pmu.ibranch)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.ibranch));
                 ];
jC11V = meas.pmu.onbus(meas.pmu.vnode);
jC11Inj = colInj;
% -------------------------------------------------------------------------

% ------------------------- Values of elements ----------------------------
vC11IijO = [     powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchOpp));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchOpp)); ];
vC11Iij =  [     powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)); 
    ];
vC11V = ones(meas.num.pV, 1);
vC11Inj = nonzeros(powsys.ybus.y(meas.pmu.loc(meas.pmu.inj), :));
% -------------------------------------------------------------------------

C11 = sparse(...
             [ iC11IijO, iC11Iij, iC11V, iC11Inj  ],... 
             [ jC11Iij0; jC11Iij; jC11V; jC11Inj  ], ...
             [ vC11IijO; vC11Iij; vC11V; vC11Inj  ] ...
         );
% -------------------------------------------------------------------------

% ------------------------- Construct C12 Matrix -------------------------- 
C12 = sparse(meas.num.pmu, powsys.num.bus);
% -------------------------------------------------------------------------
% -------------------------- Construct D Matrix ---------------------------
D = sparse(meas.num.scada, powsys.num.bus);
% -------------------------------------------------------------------------
while iter < sesettings.maxNumberOfIter  
    % --------------------------- h <-> PMU -------------------------------
    c = [   sum([ powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchOpp)), ...
            powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchOpp))] .* ...
            [ x(powsys.branch.i(-meas.pmu.loc(meas.pmu.ibranchOpp))), ...
              x(powsys.branch.j(-meas.pmu.loc(meas.pmu.ibranchOpp))) ], 2);
            sum([ powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)), ...
            powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch))] .* ...
              [ x(powsys.branch.i(meas.pmu.loc(meas.pmu.ibranch))), ...
                x(powsys.branch.j(meas.pmu.loc(meas.pmu.ibranch)))], 2);
            x(meas.pmu.onbus(meas.pmu.vnode));
            powsys.ybus.y(meas.pmu.onbus(meas.pmu.inj), :) * x(powsys.bus.busnew)
          ];
    % ---------------------------------------------------------------------
    
    % ----------------------------- d <-> SCADA ---------------------------
    d = [];
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

