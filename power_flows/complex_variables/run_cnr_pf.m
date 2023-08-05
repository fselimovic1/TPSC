function results = run_cnr_pf(solver, data)
k = 1;
converged = 0;

tic
% Initialization
if strcmp(solver.start, 'flat')
    x = ones(2 * data.nBuses, 1);
else
    x = [ data.bus(:, 8) .* exp(1i * data.bus(:, 9)); ...
          data.bus(:, 8) .* exp(1i * -data.bus(:, 9))  ];
end
  
% Compute admittance bus matrix
data.powerSystemAC = admittance_matrix(data);

% Injected current in kth iteration
Ik = complex(data.nBuses, 1);
% Bus complex power injection
Si = - (data.bus(:, 3) + 1i * data.bus(:, 4));
genbuses = data.generator(:, 1);
Si(genbuses) = Si(genbuses) + (data.generator(:, 2) + 1i * data.generator(:, 3));
% Regulated voltages
Vr = data.bus(:, 8);
Vr(genbuses) = data.generator(:, 6);

nNonZero = 0;
nInc = zeros(data.nBuses, 1);
for i = 1:data.nBuses
    if data.bus(i, 2) == 1
        nInc(i) = nnz(data.powerSystemAC.nodalMatrix(i, :));
        nNonZero = nNonZero  + 2 * (nInc(i) + 1);
    elseif data.bus(i, 2) == 2
        nInc(i) = nnz(data.powerSystemAC.nodalMatrix(i, :));
        nNonZero = nNonZero + 2 * nInc(i) + 2;
    else
        nNonZero = nNonZero + 4;
    end
end

% memory allocation
rows = zeros(nNonZero, 1);
cols = zeros(nNonZero, 1);
vals = zeros(nNonZero, 1);
g = zeros(2 * data.nBuses, 1);

cnt = 1;
for i = 1:data.nBuses
    is = [ 2 * i - 1, 2 * i ];
    if data.bus(i, 2) == 1
        Ik(i) = sum(data.powerSystemAC.nodalMatrixTranspose(:, i)...
             .* x(1:data.nBuses));
        % calculate g
        g(is) = [ x(i) * conj(Ik(i)) - Si(i), conj(x(i) * conj(Ik(i)) - Si(i)) ];
        % set indices
        nNew = nInc(i) + 1;
        idx = find(data.powerSystemAC.nodalMatrix(i, :));
        rows(cnt:cnt + 2 * nNew - 1) = [ repmat(is(1),1, nNew), ...
                                          repmat(is(2), 1, nNew) ];
        cols(cnt:cnt + 2 * nNew - 1) = [i , data.nBuses + idx, i + data.nBuses, idx ];
        cnt = cnt + 2 * nNew;
    elseif data.bus(i, 2) == 2
        Ik(i) = sum(data.powerSystemAC.nodalMatrixTranspose(:, i)...
              .* x(1:data.nBuses));
        % calculate g
        g(is) = [ real(x(i) * conj(Ik(i))) -  real(Si(i)), ...
                  x(i) * x(i + data.nBuses) - Vr(i)^2 ];
        % set indices
        nNew =  2 * nInc(i);
        rows(cnt:cnt + nNew + 1) = [ repmat(is(1), 1, nNew), ...
                                          repmat(is(2), 1, 2) ];
        idx = find(data.powerSystemAC.nodalMatrix(i, :));
        cols(cnt:cnt + nNew + 1) = [idx , data.nBuses + idx, i, data.nBuses + i ];
        cnt = cnt + nNew + 2;
    else
         % calculate g
         g(is) = [ real(x(i)) - Vr(i), imag(x(i)) ];
         % set indices
         rows(cnt:cnt + 3) = [ is(1), is(1), is(2), is(2) ];
         cols(cnt:cnt + 3) = [ i, i + data.nBuses, i , i + data.nBuses ];
         cnt = cnt + 4;
    end
end
J = sparse(rows, cols, vals, 2 * data.nBuses, 2 * data.nBuses);

while k < solver.maxNumberOfIter
    % compute nonzeros of J
    cnt = 1;
    for i = 1:data.nBuses
        if data.bus(i, 2) == 1
            nNew = nInc(i) + 1;
            vals(cnt:cnt + 2 * nNew - 1) = [ conj(Ik(i)); nonzeros(x(i) .* ...
                        conj(data.powerSystemAC.nodalMatrixTranspose(:, i))); ...
                        Ik(i); nonzeros(x(i + data.nBuses) .* ...
                        data.powerSystemAC.nodalMatrixTranspose(:, i))];
            cnt = cnt + 2 * nNew;
        elseif data.bus(i, 2) == 2
            nNew =  2 * nInc(i);
            idx = find(data.powerSystemAC.nodalMatrix(i, :));
            iIdx = find(idx == i);
            valsS = [ zeros(iIdx - 1, 1); conj(Ik(i)); zeros(nInc(i) - iIdx, 1);...
                      nonzeros(x(i) .* conj(data.powerSystemAC.nodalMatrixTranspose(:, i))) ];
            vals(cnt:cnt +  nNew + 1) = [ 1/2 .* (valsS + conj([valsS(nInc(i) + 1:2*nInc(i))...
                                          ;valsS(1:nInc(i))]));  ...
                                          x(i + data.nBuses); x(i) ];
            cnt = cnt + nNew + 2;
        else 
            vals(cnt:cnt + 3) = [ 1/2, 1/2, -1i/2, 1i/2 ];
            cnt = cnt + 4;
        end
    end
    for i = 1:nNonZero
        J(rows(i), cols(i)) = vals(i);
    end
    % calculate new iterates
    x = x - J \ g;
    
    % compute residuals
    for i = 1:data.nBuses
        if data.bus(i, 2) == 1
            Ik(i) = sum(data.powerSystemAC.nodalMatrixTranspose(:, i)...
             .* x(1:data.nBuses));
            g([2 * i - 1, 2 * i]) = [ x(i) * conj(Ik(i)) - Si(i), ...
                                     conj(x(i) * conj(Ik(i)) - Si(i)) ];
        elseif data.bus(i, 2) == 2
            Ik(i) = sum(data.powerSystemAC.nodalMatrixTranspose(:, i)...
             .* x(1:data.nBuses));
            g([2 * i - 1, 2 * i]) = [ real(x(i) * conj(Ik(i)))...
                                      - real(Si(i)), x(i) * x(i + data.nBuses) - Vr(i)^2 ];
        else
            g([2 * i - 1, 2 * i]) = [ real(x(i)) - Vr(i), imag(x(i)) ];
        end
    end
    
    % check convergence
    if max(abs(g)) < solver.eps
        converged = 1;
        break;
    end
    k = k + 1;
end
algtime = toc;
% post-processing
if solver.postprocess
    tic
    results = postprocess_acpf(data.powerSystemAC, x(1:data.nBuses));
    results.ppt = toc;
end
results.at = algtime;
results.iter = k;
results.converged = converged;
results.sys = data.case;
results.method = 'Newton-Raphson in Complex Variables';
end

