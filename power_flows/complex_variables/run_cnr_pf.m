function [ V, iter, converged, method ] = run_cnr_pf(pfsettings, powersystem, ybus, colsYbus, num)
iter = 1;
converged = 0;

% Initialization
if pfsettings.flatStart
    x = ones(2 * num.bus, 1);
else
    x = [ powersystem.Vmi .* exp(1i * powersystem.Vai); ...
          powersystem.Vmi .* exp(1i * -powersystem.Vai)  ];
end

% Bus complex power injection
Si = -(powersystem.Pload + 1i * powersystem.Qload);
Si(powersystem.genbuses) = Si(powersystem.genbuses) + (powersystem.Pgeni + 1i * powersystem.Qgeni);

% Regulated voltages
Vr = powersystem.Vmi;
Vr(powersystem.genbuses) = powersystem.Vgen;
Vsl = powersystem.Vmi(num.islack) * exp(1i * powersystem.Vai(num.islack));

% important numbers
nInc = sum(ybus.nodalMatrix ~= 0, 2);
nNonZero = 4 + 2 * sum(nInc(powersystem.isPQ | powersystem.isPV) + 1);
nPQ = sum(powersystem.isPQ);
nPV = sum(powersystem.isPV);
pq = powersystem.bus(powersystem.isPQ);
pv = powersystem.bus(powersystem.isPV);

% memory allocation
rows = zeros(nNonZero, 1);
cols = zeros(nNonZero, 1);
vals = zeros(nNonZero, 1);

cnt = 1;
% PQ buses
for i = 1:nPQ
    pqbus = pq(i);
    % set indices
    nNew = nInc(pqbus) + 1;
    rows(cnt:cnt + 2 * nNew - 1) = [ repmat(i, 1, nNew), ...
                                     repmat(nPQ +  i, 1, nNew) ];
    cols(cnt:cnt + 2 * nNew - 1) = [ pqbus; num.bus + colsYbus{pqbus}; pqbus + num.bus; colsYbus{pqbus} ];
    cnt = cnt + 2 * nNew;
end
jpq = cnt - 1;
% PV buses
for i = 1:nPV
    pvbus = pv(i);
    % set indices
    nNew =  2 * nInc(pvbus);
    rows(cnt:cnt + nNew + 1) = [ repmat(2 * nPQ + i, 1, nNew), ...
                                 repmat(2 * nPQ + nPV + i, 1, 2) ];
    cols(cnt:cnt + nNew + 1) = [colsYbus{pvbus}; num.bus + colsYbus{pvbus}; pvbus; num.bus + pvbus ];
    cnt = cnt + nNew + 2;
end
auxPV = zeros(nPV, num.bus);
linpvidx = (pv - 1) * nPV + (1:nPV);
jpv = cnt - jpq - 1;
% SLACK bus
rows(cnt:cnt + 3) = [ 2 * num.bus - 1, 2 * num.bus - 1, 2 * num.bus, 2 * num.bus ];
cols(cnt:cnt + 3) = [ num.islack, num.islack + num.bus, num.islack , num.islack + num.bus ];

J = sparse(rows, cols, zeros(nNonZero, 1), 2 * num.bus, 2 * num.bus);
linidx = (cols - 1) * 2 * num.bus + rows;
g = zeros(2 * num.bus, 1);

Ik = ybus.nodalMatrix * x(powersystem.bus) + eps; % + epse -> to avoid zeros
% function values at the current state
g(1:2 * nPQ) = [ x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ),...
    conj(x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ))];
g(2 * nPQ + 1:2 * nPQ + 2 * nPV) = [ real(x(powersystem.isPV) .* ...
    conj(Ik(powersystem.isPV)) - Si(powersystem.isPV)), x(powersystem.isPV)...
    .* x([ false(num.bus, 1), powersystem.isPV ]) - Vr(powersystem.isPV).^2 ];
g(2 * nPQ + 2 * nPV + 1:end) = [ real(x(num.islack)) - real(Vsl),...
    imag(x(num.islack)) - imag(Vsl) ];

while iter < pfsettings.maxNumberOfIter
	% compute nonzeros of J
	% PQ buses
    vals(1:jpq/2) = nonzeros(reshape([conj(Ik(powersystem.isPQ)), x(powersystem.bus).' .* ...
	conj(ybus.nodalMatrix(powersystem.isPQ, :))], [], 1));
    vals(jpq/2 + 1:jpq) = nonzeros(reshape([x(powersystem.bus)' .* ...
	ybus.nodalMatrix(powersystem.isPQ, :), Ik(powersystem.isPQ)], [], 1));
    % PV buses
	auxPV(linpvidx) = Ik(pv); 
    vals(jpq + 1:jpq + jpv - 2 * nPV) = nonzeros(reshape(1/2 .* ([conj(auxPV), x(powersystem.bus).' .* ...
	conj(ybus.nodalMatrix(powersystem.isPV, :))] + [ x(powersystem.bus)' .* ...
	ybus.nodalMatrix(powersystem.isPV, :), auxPV]), [], 1));
    vals(jpq + jpv - 2 * nPV + 1:jpq + jpv) = reshape([ conj(x(powersystem.bus(powersystem.isPV))), x(powersystem.bus(powersystem.isPV)) ], [], 1);
	% SLACK bus
    vals(jpq + jpv + 1:end) =  [ 1/2, 1/2, -1i/2, 1i/2 ];
    
    J(linidx) = vals;
    % calculate new iterates
    x = x - J \ g;
    
   Ik = ybus.nodalMatrix * x(powersystem.bus) + eps; % + epse -> to avoid zeros
    % function values at the current state
    g(1:2 * nPQ) = [ x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ),...
        conj(x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ))];
    g(2 * nPQ +1:2 * nPQ + 2 * nPV) = [ real(x(powersystem.isPV) .* ...
        conj(Ik(powersystem.isPV)) - Si(powersystem.isPV)), x(powersystem.isPV)...
        .* x([ false(num.bus, 1), powersystem.isPV ]) - Vr(powersystem.isPV).^2 ];
    g(2 * nPQ + 2 * nPV + 1:end) = [ real(x(num.islack)) - real(Vsl),...
        imag(x(num.islack)) - imag(Vsl) ];
    
    % check convergence
    if max(abs(g)) < pfsettings.eps
        converged = 1;
        break;
    end
    iter = iter + 1;
end
V = x(powersystem.bus);
method = 'Newton-Raphson in Complex Variables';
end

