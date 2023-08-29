function [ V, iter, converged, method ] = run_cnr_pf(powsys, pfsettings)
nPQ = numel(powsys.bus.pq);
nPV = numel(powsys.bus.pv);
iter = 1;
converged = 0;

% ---------------------- Initialize state variables -----------------------
if pfsettings.flatStart
    x = ones(2 * powsys.num.bus, 1);
else
    x = [ powsys.bus.Vmi .* exp(1i * powsys.bus.Vai); ...
          powsys.bus.Vmi .* exp(1i * -powsys.bus.Vai)  ];
end
% -------------------------------------------------------------------------

% ----------------------- Complex power injection -------------------------
Si = -(powsys.bus.Pload + 1i * powsys.bus.Qload) + powsys.bus.Pgen + 1i * powsys.bus.Qgen;
% -------------------------------------------------------------------------

% --------------------- Compute function g values -------------------------
Ik = powsys.ybus.y * x(powsys.bus.busnew);
S = x(powsys.bus.busnew) .* conj(Ik);
g = [ S(powsys.bus.pq) -  Si(powsys.bus.pq); ...
      conj(S(powsys.bus.pq) -  Si(powsys.bus.pq)); ...
      real(S(powsys.bus.pv)) - real(Si(powsys.bus.pv)); ...
      x(powsys.bus.pv) .* conj(x(powsys.bus.pv)) - powsys.bus.Vmi(powsys.bus.pv) .^ 2; ...
      real(x(powsys.num.islack)) - real(powsys.bus.Vmi(powsys.num.islack) ...
                            * exp(1i * powsys.bus.Vai(powsys.num.islack))); ...
      imag(x(powsys.num.islack)) - imag(powsys.bus.Vmi(powsys.num.islack) ...
                            * exp(1i * powsys.bus.Vai(powsys.num.islack)));
       ];
% -------------------------------------------------------------------------

% ----------------------- J matrix indices --------------------------------
[ c12i, c12j ] = find(powsys.ybus.y(powsys.bus.pq, :));
[ dyi, dyj ] = find(powsys.ybus.yij(powsys.bus.pv, :));
di = [ (1:nPV)'; dyi; nPV + (1:nPV)';  2 * nPV + 1; 2 * nPV + 2 ];
dj = [ powsys.bus.pv; dyj; powsys.bus.pv;  powsys.num.islack; powsys.num.islack ];
% -------------------------------------------------------------------------

while iter < pfsettings.maxNumberOfIter
    % ------------------- Calculate J matrix entries ----------------------
    % ------------------------ MATRIX C11 ---------------------------------
    C11 = sparse(1:nPQ, powsys.bus.pq, conj(Ik(powsys.bus.pq)), nPQ, powsys.num.bus);
    % ---------------------------------------------------------------------
    % ------------------------ MATRIX C12 ---------------------------------
	C12 = sparse(c12i, c12j, nonzeros(x(powsys.bus.busnew).' .* conj(powsys.ybus.y(powsys.bus.pq, :))), nPQ, powsys.num.bus);
    % ---------------------------------------------------------------------
    % ------------------------- MATRIX D ----------------------------------
    D = sparse(di, dj, [ 1/2 .* [ conj(Ik(powsys.bus.pv)) + x(powsys.bus.pv + powsys.num.bus) .*  ...
        powsys.ybus.ydiag(powsys.bus.pv); ...
        nonzeros(x(powsys.bus.busnew)' .* powsys.ybus.yij( powsys.bus.pv, :)) ]; ...
        conj(x(powsys.bus.pv)); 1/2; -1i/2 ], 2 * nPV + 2, powsys.num.bus);
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    J = [ C11,        C12;
          conj(C12)  conj(C11);
            D          conj(D) ];
    % ---------------------- Calculate new iterates -----------------------
    x = x - J \ g;
    % ---------------------------------------------------------------------
    
    % --------------------- Compute function g values -------------------------
    Ik = powsys.ybus.y * x(powsys.bus.busnew);
    S = x(powsys.bus.busnew) .* conj(Ik);
    g = [ S(powsys.bus.pq) -  Si(powsys.bus.pq); ...
          conj(S(powsys.bus.pq) -  Si(powsys.bus.pq)); ...
          real(S(powsys.bus.pv)) - real(Si(powsys.bus.pv)); ...
          x(powsys.bus.pv) .* conj(x(powsys.bus.pv)) - powsys.bus.Vmi(powsys.bus.pv) .^ 2; ...
          real(x(powsys.num.islack)) - real(powsys.bus.Vmi(powsys.num.islack) ...
                                * exp(1i * powsys.bus.Vai(powsys.num.islack))); ...
          imag(x(powsys.num.islack)) - imag(powsys.bus.Vmi(powsys.num.islack) ...
                                * exp(1i * powsys.bus.Vai(powsys.num.islack)));
           ];
    % -------------------------------------------------------------------------

    %  -------------------- Convergence check -----------------------------
    if max(abs(g)) < pfsettings.eps
        converged = 1;
        break;
    end
    iter = iter + 1;
    % ---------------------------------------------------------------------
end
V = x(powsys.bus.busnew);
method = 'Newton-Raphson in Complex Variables';
end

