function [ V, iter, converged, method ] = run_cnr_pf(pfsettings, powersystem, ybus, num)
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

nInc = sum(ybus.nodalMatrix ~= 0, 2);
nNonZero = 4 + 2 * sum(nInc(powersystem.isPQ | powersystem.isPV) + 1);

% memory allocation
rows = zeros(nNonZero, 1);
cols = zeros(nNonZero, 1);
g = zeros(2 * num.bus, 1);

cnt = 1;

Ik = ybus.nodalMatrix * x(powersystem.bus);
% function values at the current state
g(repmat(powersystem.isPQ, 2, 1)) = [ x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ),...
    conj(x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ))];
g(repmat(powersystem.isPV, 2, 1)) = [ real(x(powersystem.isPV) .* ...
    conj(Ik(powersystem.isPV)) - Si(powersystem.isPV)), x(powersystem.isPV)...
    .* x([ false(num.bus, 1), powersystem.isPV ]) - Vr(powersystem.isPV).^2 ];
g([num.islack, num.islack + num.bus]) = [ real(x(num.islack)) - real(Vsl),...
    imag(x(num.islack)) - imag(Vsl) ];
for i = 1:num.bus
    if powersystem.bustype(i) == 1
        % set indices
        nNew = nInc(i) + 1;
        idx = find(ybus.nodalMatrix(i, :));
        rows(cnt:cnt + 2 * nNew - 1) = [ repmat(i, 1, nNew), ...
                                          repmat(num.bus +  i, 1, nNew) ];
        cols(cnt:cnt + 2 * nNew - 1) = [i , num.bus + idx, i + num.bus, idx ];
        cnt = cnt + 2 * nNew;
    elseif powersystem.bustype(i) == 2
        % set indices
        nNew =  2 * nInc(i);
        rows(cnt:cnt + nNew + 1) = [ repmat(i, 1, nNew), ...
                                     repmat(num.bus + i, 1, 2) ];
        idx = find(ybus.nodalMatrix(i, :));
        cols(cnt:cnt + nNew + 1) = [idx , num.bus + idx, i, num.bus + i ];
        cnt = cnt + nNew + 2;
    else
         % set indices
         rows(cnt:cnt + 3) = [ i, i, i + num.bus, i + num.bus ];
         cols(cnt:cnt + 3) = [ i, i + num.bus, i , i + num.bus ];
         cnt = cnt + 4;
    end
end
J = sparse(rows, cols, zeros(nNonZero, 1), 2 * num.bus, 2 * num.bus);

while iter < pfsettings.maxNumberOfIter
    % compute nonzeros of J
    cnt = 1;
    for i = 1:num.bus
        if powersystem.bustype(i) == 1
            nNew = nInc(i) + 1;
            J(i, cols(cnt:cnt + nNew - 1)) = [ conj(Ik(i)); nonzeros(x(i) .* ...
                        conj(ybus.nodalMatrixTranspose(:, i))) ];
            cnt = cnt + nNew;
            J(num.bus + i, cols(cnt:cnt + nNew - 1)) = [ Ik(i); nonzeros(x(i + num.bus) .* ...
                                     ybus.nodalMatrixTranspose(:, i)) ];
            cnt = cnt + nNew;
        elseif powersystem.bustype(i) == 2
            nNew =  2 * nInc(i);
            idx = find(ybus.nodalMatrix(i, :));
            iIdx = find(idx == i);
            valsS = [ zeros(iIdx - 1, 1); conj(Ik(i)); zeros(nInc(i) - iIdx, 1);...
                      nonzeros(x(i) .* conj(ybus.nodalMatrixTranspose(:, i))) ];
            J(i, cols(cnt:cnt +  nNew - 1)) = [ 1/2 .* (valsS + conj([valsS(nInc(i) + 1:2 * nInc(i))...
                                          ; valsS(1:nInc(i))])) ];
            cnt = cnt + nNew;
            J(num.bus + i, [i, num.bus + i]) = [ x(i + num.bus); x(i) ];                       
            cnt = cnt + 2;
        else 
            J(i, [i, num.bus + i]) = [ 1/2, 1/2 ];
            J(num.bus + i, [i, num.bus + i]) = [ -1i/2, 1i/2 ];
            cnt = cnt + 4;
        end
    end
    
    % calculate new iterates
    x = x - J \ g;
    
    Ik = ybus.nodalMatrix * x(powersystem.bus);
    % function values at the current state
    g(repmat(powersystem.isPQ, 2, 1)) = [ x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ),...
        conj(x(powersystem.isPQ) .* conj(Ik(powersystem.isPQ)) - Si(powersystem.isPQ))];
    g(repmat(powersystem.isPV, 2, 1)) = [ real(x(powersystem.isPV) .* ...
        conj(Ik(powersystem.isPV)) - Si(powersystem.isPV)), x(powersystem.isPV)...
        .* x([ false(num.bus, 1), powersystem.isPV ]) - Vr(powersystem.isPV).^2 ];
    g([num.islack, num.islack + num.bus]) = [ real(x(num.islack)) - real(Vsl),...
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

