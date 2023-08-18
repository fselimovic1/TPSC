function results = run_cnr_pf(pfsettings, data, varargin)
trackingmode = 0;
if nargin - 2 > 0
    trackingmode = 1;
    dynsettings = varargin{1};
end
k = 1;
converged = 0;

% renumber data if needed
bus = 1:data.nBuses;
busi = data.bus(:, 1);
genbuses = data.generator(:, 1);
branchi = data.branch(:, 1);
branchj = data.branch(:, 2);
toRenumber = any(data.bus(:, 1) ~= bus);
if toRenumber
    genbuses = renumbering(genbuses, busi, bus);
    branchi = renumbering(branchi, busi, bus);
    branchj = renumbering(branchj, busi, bus);
end


tic
% Initialization
if strcmp(pfsettings.start, 'flat')
    x = ones(2 * data.nBuses, 1);
else
    x = [ data.bus(:, 8) .* exp(1i * data.bus(:, 9)); ...
          data.bus(:, 8) .* exp(1i * -data.bus(:, 9))  ];
end
  
% Compute admittance bus matrix
if trackingmode && dynsettings.f ~= -data.fn
    data.powerSystemAC = admittance_matrix(data, branchi, branchj, dynsettings);
else
    data.powerSystemAC = admittance_matrix(data, branchi, branchj);
end

% Injected current in kth iteration
Ik = complex(data.nBuses, 1);

% Bus complex power injection
Si = - (data.bus(:, 3) + 1i * data.bus(:, 4));
if trackingmode && dynsettings.loadNo 
    if dynsettings.loadNo == -1
        Si = Si .* dynsettings.load;
    else
        Si(dynsettings.loadNo) = Si(dynsettings.loadNo) + ...
                                 complex(dynsettings.load(1), dynsettings.load(2));  
    end
end
Si(genbuses) = Si(genbuses) + (data.generator(:, 2) + 1i * data.generator(:, 3));

% Regulated voltages
Vr = data.bus(:, 8);
Vr(genbuses) = data.generator(:, 6);
Vsl = data.bus(data.slackNo, 8) * exp(1i * data.bus(data.slackNo, 9));

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
g = zeros(2 * data.nBuses, 1);

cnt = 1;
for i = 1:data.nBuses
    if data.bus(i, 2) == 1
        Ik(i) = sum(data.powerSystemAC.nodalMatrixTranspose(:, i)...
             .* x(1:data.nBuses));
        % calculate g
        g([2 * i - 1, 2 * i]) = [ x(i) * conj(Ik(i)) - Si(i), conj(x(i) * conj(Ik(i)) - Si(i)) ];
        % set indices
        nNew = nInc(i) + 1;
        idx = find(data.powerSystemAC.nodalMatrix(i, :));
        rows(cnt:cnt + 2 * nNew - 1) = [ repmat(2 * i - 1,1, nNew), ...
                                          repmat(2 * i, 1, nNew) ];
        cols(cnt:cnt + 2 * nNew - 1) = [i , data.nBuses + idx, i + data.nBuses, idx ];
        cnt = cnt + 2 * nNew;
    elseif data.bus(i, 2) == 2
        Ik(i) = sum(data.powerSystemAC.nodalMatrixTranspose(:, i)...
              .* x(1:data.nBuses));
        % calculate g
        g([2 * i - 1, 2 * i]) = [ real(x(i) * conj(Ik(i))) -  real(Si(i)), ...
                  x(i) * x(i + data.nBuses) - Vr(i)^2 ];
        % set indices
        nNew =  2 * nInc(i);
        rows(cnt:cnt + nNew + 1) = [ repmat(2 * i - 1, 1, nNew), ...
                                     repmat(2 * i, 1, 2) ];
        idx = find(data.powerSystemAC.nodalMatrix(i, :));
        cols(cnt:cnt + nNew + 1) = [idx , data.nBuses + idx, i, data.nBuses + i ];
        cnt = cnt + nNew + 2;
    else
         % calculate g
         g([2 * i - 1, 2 * i]) = [ real(x(i)) - real(Vsl), imag(x(i)) - imag(Vsl) ];
         % set indices
         rows(cnt:cnt + 3) = [ 2 * i - 1, 2 * i - 1, 2 * i, 2 * i ];
         cols(cnt:cnt + 3) = [ i, i + data.nBuses, i , i + data.nBuses ];
         cnt = cnt + 4;
    end
end
J = sparse(rows, cols, zeros(nNonZero, 1), 2 * data.nBuses, 2 * data.nBuses);

while k < pfsettings.maxNumberOfIter
    % compute nonzeros of J
    cnt = 1;
    for i = 1:data.nBuses
        if data.bus(i, 2) == 1
            nNew = nInc(i) + 1;
            J(2 * i - 1, cols(cnt:cnt + nNew - 1)) = [ conj(Ik(i)); nonzeros(x(i) .* ...
                        conj(data.powerSystemAC.nodalMatrixTranspose(:, i))) ];
            cnt = cnt + nNew;
            J(2 * i, cols(cnt:cnt + nNew - 1)) = [ Ik(i); nonzeros(x(i + data.nBuses) .* ...
                                     data.powerSystemAC.nodalMatrixTranspose(:, i)) ];
            cnt = cnt + nNew;
        elseif data.bus(i, 2) == 2
            nNew =  2 * nInc(i);
            idx = find(data.powerSystemAC.nodalMatrix(i, :));
            iIdx = find(idx == i);
            valsS = [ zeros(iIdx - 1, 1); conj(Ik(i)); zeros(nInc(i) - iIdx, 1);...
                      nonzeros(x(i) .* conj(data.powerSystemAC.nodalMatrixTranspose(:, i))) ];
            J(2 * i - 1, cols(cnt:cnt +  nNew - 1)) = [ 1/2 .* (valsS + conj([valsS(nInc(i) + 1:2 * nInc(i))...
                                          ; valsS(1:nInc(i))])) ];
            cnt = cnt + nNew;
            J(2 * i, [i, data.nBuses + i]) = [ x(i + data.nBuses); x(i) ];                       
            cnt = cnt + 2;
        else 
            J(2 * i - 1, [i, data.nBuses + i]) = [ 1/2, 1/2 ];
            J(2 * i, [i, data.nBuses + i]) = [ -1i/2, 1i/2 ];
            cnt = cnt + 4;
        end
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
            g([2 * i - 1, 2 * i]) = [ real(x(i)) - real(Vsl), imag(x(i)) - imag(Vsl) ];
        end
    end
    
    % check convergence
    if max(abs(g)) < pfsettings.eps
        converged = 1;
        break;
    end
    k = k + 1;
end
algtime = toc;
% post-processing
if pfsettings.postprocess
    tic
    if trackingmode
        v = abs(x(1:data.nBuses)) * exp(1i * (angle(x(data.nBuses)) + dynsettings.thetaslack));
    else
        v = x(1:data.nBuses);
    end
    results = postprocess_acpf(data.powerSystemAC, branchi, branchj, v);
    results.Pload = data.bus(:, 3);
    results.Qload = data.bus(:, 4);
    results.Pgen(genbuses) = data.generator(:, 2);
    results.Qgen(genbuses) = results.Qi(genbuses)  + results.Qload(genbuses);
    results.ppt = toc;
end
results.at = algtime;
results.iter = k;
results.converged = converged;
results.sys = data.case;
results.method = 'Newton-Raphson in Complex Variables';
end

