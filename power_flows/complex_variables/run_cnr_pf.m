function results = run_cnr_pf(pfsettings, data, varargin)
trackingmode = 0;
if nargin - 2 > 0
    trackingmode = 1;
    dynsettings = varargin{1};
end
k = 1;
converged = 0;

% number of:
num.bus = size(data.bus, 1);
num.branch = size(data.branch, 1);
num.gen = size(data.generator, 1);

% read the power system data
bus = 1:num.bus;
busi = data.bus(:, 1);
branchi = data.branch(:, 1);
branchj = data.branch(:, 2);
resistance = data.branch(:, 3);
reactance = data.branch(:, 4);
charging = data.branch(:, 5);
transturns = data.branch(:, 9);
transphase = data.branch(:, 10);
branchstatus = data.branch(:, 11) == 1;
bustype = data.bus(:, 2);
num.islack= find(bustype == 3);
Pload = data.bus(:, 3);
Qload = data.bus(:, 4);
Gshunt = data.bus(:, 5);
Bshunt = data.bus(:, 6);
Vmi = data.bus(:, 8);
Vai = data.bus(:, 9);
genbuses = data.generator(:, 1);
Pgeni = data.generator(:, 2);
Qgeni = data.generator(:, 3);
% Qgmax, Qgmin
Vgen = data.generator(:, 6);

% renumber data if needed
toRenumber = any(data.bus(:, 1) ~= bus);
if toRenumber
    genbuses = renumbering(genbuses, busi, bus);
    branchi = renumbering(branchi, busi, bus);
    branchj = renumbering(branchj, busi, bus);
end

tic
% Initialization
if pfsettings.flatStart
    x = ones(2 * num.bus, 1);
else
    x = [ Vmi .* exp(1i * Vai); ...
          Vmi .* exp(1i * -Vai)  ];
end
  
% Compute admittance bus matrix
powerSystemAC = admittance_matrix(num, bus, branchi, branchj, resistance, reactance, charging, transturns, transphase, branchstatus, Gshunt, Bshunt);

% Injected current in kth iteration
Ik = complex(num.bus, 1);

% Bus complex power injection
Si = - (Pload + 1i * Qload);
if trackingmode && dynsettings.loadNo 
    if dynsettings.loadNo == -1
        Si = Si .* dynsettings.load;
    else
        Si(dynsettings.loadNo) = Si(dynsettings.loadNo) + ...
                                 complex(dynsettings.load(1), dynsettings.load(2));  
    end
end
Si(genbuses) = Si(genbuses) + (Pgeni + 1i * Qgeni);

% Regulated voltages
Vr = Vmi;
Vr(genbuses) = Vgen;
Vsl = Vmi(num.islack) * exp(1i * Vai(num.islack));

nNonZero = 0;
nInc = zeros(num.bus, 1);
for i = 1:num.bus
    if bustype(i) == 1
        nInc(i) = nnz(powerSystemAC.nodalMatrix(i, :));
        nNonZero = nNonZero  + 2 * (nInc(i) + 1);
    elseif bustype(i) == 2
        nInc(i) = nnz(powerSystemAC.nodalMatrix(i, :));
        nNonZero = nNonZero + 2 * nInc(i) + 2;
    else
        nNonZero = nNonZero + 4;
    end
end

% memory allocation
rows = zeros(nNonZero, 1);
cols = zeros(nNonZero, 1);
g = zeros(2 * num.bus, 1);

cnt = 1;
for i = 1:num.bus
    if bustype(i) == 1
        Ik(i) = sum(powerSystemAC.nodalMatrixTranspose(:, i)...
             .* x(bus));
        % calculate g
        g([2 * i - 1, 2 * i]) = [ x(i) * conj(Ik(i)) - Si(i), conj(x(i) * conj(Ik(i)) - Si(i)) ];
        % set indices
        nNew = nInc(i) + 1;
        idx = find(powerSystemAC.nodalMatrix(i, :));
        rows(cnt:cnt + 2 * nNew - 1) = [ repmat(2 * i - 1,1, nNew), ...
                                          repmat(2 * i, 1, nNew) ];
        cols(cnt:cnt + 2 * nNew - 1) = [i , num.bus + idx, i + num.bus, idx ];
        cnt = cnt + 2 * nNew;
    elseif bustype(i) == 2
        Ik(i) = sum(powerSystemAC.nodalMatrixTranspose(:, i)...
              .* x(bus));
        % calculate g
        g([2 * i - 1, 2 * i]) = [ real(x(i) * conj(Ik(i))) -  real(Si(i)), ...
                  x(i) * x(i + num.bus) - Vr(i)^2 ];
        % set indices
        nNew =  2 * nInc(i);
        rows(cnt:cnt + nNew + 1) = [ repmat(2 * i - 1, 1, nNew), ...
                                     repmat(2 * i, 1, 2) ];
        idx = find(powerSystemAC.nodalMatrix(i, :));
        cols(cnt:cnt + nNew + 1) = [idx , num.bus + idx, i, num.bus + i ];
        cnt = cnt + nNew + 2;
    else
         % calculate g
         g([2 * i - 1, 2 * i]) = [ real(x(i)) - real(Vsl), imag(x(i)) - imag(Vsl) ];
         % set indices
         rows(cnt:cnt + 3) = [ 2 * i - 1, 2 * i - 1, 2 * i, 2 * i ];
         cols(cnt:cnt + 3) = [ i, i + num.bus, i , i + num.bus ];
         cnt = cnt + 4;
    end
end
J = sparse(rows, cols, zeros(nNonZero, 1), 2 * num.bus, 2 * num.bus);

while k < pfsettings.maxNumberOfIter
    % compute nonzeros of J
    cnt = 1;
    for i = 1:num.bus
        if bustype(i) == 1
            nNew = nInc(i) + 1;
            J(2 * i - 1, cols(cnt:cnt + nNew - 1)) = [ conj(Ik(i)); nonzeros(x(i) .* ...
                        conj(powerSystemAC.nodalMatrixTranspose(:, i))) ];
            cnt = cnt + nNew;
            J(2 * i, cols(cnt:cnt + nNew - 1)) = [ Ik(i); nonzeros(x(i + num.bus) .* ...
                                     powerSystemAC.nodalMatrixTranspose(:, i)) ];
            cnt = cnt + nNew;
        elseif bustype(i) == 2
            nNew =  2 * nInc(i);
            idx = find(powerSystemAC.nodalMatrix(i, :));
            iIdx = find(idx == i);
            valsS = [ zeros(iIdx - 1, 1); conj(Ik(i)); zeros(nInc(i) - iIdx, 1);...
                      nonzeros(x(i) .* conj(powerSystemAC.nodalMatrixTranspose(:, i))) ];
            J(2 * i - 1, cols(cnt:cnt +  nNew - 1)) = [ 1/2 .* (valsS + conj([valsS(nInc(i) + 1:2 * nInc(i))...
                                          ; valsS(1:nInc(i))])) ];
            cnt = cnt + nNew;
            J(2 * i, [i, num.bus + i]) = [ x(i + num.bus); x(i) ];                       
            cnt = cnt + 2;
        else 
            J(2 * i - 1, [i, num.bus + i]) = [ 1/2, 1/2 ];
            J(2 * i, [i, num.bus + i]) = [ -1i/2, 1i/2 ];
            cnt = cnt + 4;
        end
    end
    
    % calculate new iterates
    x = x - J \ g;
    
    % compute residuals
    for i = 1:num.bus
        if bustype(i) == 1
            Ik(i) = sum(powerSystemAC.nodalMatrixTranspose(:, i)...
             .* x(bus));
            g([2 * i - 1, 2 * i]) = [ x(i) * conj(Ik(i)) - Si(i), ...
                                     conj(x(i) * conj(Ik(i)) - Si(i)) ];
        elseif bustype(i) == 2
            Ik(i) = sum(powerSystemAC.nodalMatrixTranspose(:, i)...
             .* x(bus));
            g([2 * i - 1, 2 * i]) = [ real(x(i) * conj(Ik(i)))...
                                      - real(Si(i)), x(i) * x(i + num.bus) - Vr(i)^2 ];
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
        v = abs(x(bus)) * exp(1i * (angle(x(num.bus)) + dynsettings.thetaslack));
    else
        v = x(bus);
    end
    results = postprocess_acpf(powerSystemAC, branchi, branchj, v);
    results.Pload = Pload;
    results.Qload = Qload;
    results.Pgen(genbuses) = Pgeni;
    results.Qgen(genbuses) = results.Qi(genbuses)  + results.Qload(genbuses);
    results.ppt = toc;
end
results.at = algtime;
results.iter = k;
results.converged = converged;
results.sys = data.case;
results.method = 'Newton-Raphson in Complex Variables';
end

