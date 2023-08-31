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
      meas.pmu.m(meas.pmu.ibranchOpp) .* exp(-1i .* meas.pmu.a(meas.pmu.ibranchOpp));
      meas.pmu.m(meas.pmu.ibranch) .* exp(-1i .* meas.pmu.a(meas.pmu.ibranch));
      meas.pmu.m(meas.pmu.vnode) .* exp(-1i .* meas.pmu.a(meas.pmu.vnode));
      meas.scada.m
    ];
% -------------------------------------------------------------------------

% ------------------------ Measurements' weights --------------------------
W = sparse(1:2 * meas.num.pmu + meas.num.scada, 1:2 * meas.num.pmu + meas.num.scada, ...
           [ str2double(sesettings.mweights(2)) .* ones( 2 * meas.num.pmu, 1);...
           ones(meas.num.scada, 1) ]);
% -------------------------------------------------------------------------

% ------------------------- Construct C11 matrix --------------------------
C11 = sparse(...
             [ (1:meas.num.pIijO ), (1:meas.num.pIijO ), (meas.num.pIijO  +...
               (1:meas.num.pIij)), (meas.num.pIijO  + (1:meas.num.pIij)),...
                (meas.num.pIijO + meas.num.pIij + (1:meas.num.pV))  ],... 
             [   powsys.branch.i(-meas.pmu.loc(meas.pmu.ibranchOpp)); 
                 powsys.branch.j(-meas.pmu.loc(meas.pmu.ibranchOpp)); 
                 powsys.branch.i(meas.pmu.loc(meas.pmu.ibranch)); 
                 powsys.branch.j(meas.pmu.loc(meas.pmu.ibranch));
                 meas.pmu.onbus(meas.pmu.vnode) ], ...
             [   powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchOpp));...
                 powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchOpp)); ...
                 powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch));...
                 powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch)); 
                 ones(meas.num.pV, 1)] ...
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
%     c = [ diag([ powsys.ybus.tofrom(-meas.pmu.loc(meas.pmu.ibranchOpp)), ...
%                  powsys.ybus.toto(-meas.pmu.loc(meas.pmu.ibranchOpp))] * ...
%                [ x(powsys.branch.i(-meas.pmu.loc(meas.pmu.ibranchOpp))).';
%                  x(powsys.branch.j(-meas.pmu.loc(meas.pmu.ibranchOpp))).' ] );
%          diag([ powsys.ybus.fromfrom(meas.pmu.loc(meas.pmu.ibranch)), ...
%          powsys.ybus.fromto(meas.pmu.loc(meas.pmu.ibranch))] * ...
%        [ x(powsys.branch.i(meas.pmu.loc(meas.pmu.ibranch))).';
%          x(powsys.branch.j(meas.pmu.loc(meas.pmu.ibranch))).' ] );
%          x(meas.pmu.onbus(meas.pmu.vnode))
%           ];
    % ---------------------------------------------------------------------
    
    % ----------------------------- d <-> SCADA ---------------------------
    d = [];
    % ---------------------------------------------------------------------
    H = [  C11         C12
         conj(C12)   conj(C11)
           D      conj(D)];
%     h = [c; conj(c); d];
    
    x = H \ z;
    converged = 1;
    break;
%     % Solve
%     r = z - h;
%     dx = (H' * W * H) \ (H' * W * r);
%     x = x + dx; 
%     if max(abs(dx)) < sesettings.eps
%         converged = 1;
%         break
%     else
%         iter = iter + 1;
%     end
end
% ----------------------- Estimator results -------------------------------
Vc = x(powsys.bus.busnew);
info.nonZerosInH = nnz(H);
info.redundancy = (2 * meas.num.pmu + meas.num.scada)/(2 * powsys.num.bus);
% -------------------------------------------------------------------------
end

