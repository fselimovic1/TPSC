function results = run_cnr_pf(solver, data)
iNo = 1;
nPQ = sum(data.bus(:, 2) == 1);
% Compute admittance bus matrix
data.powerSystemAC = admittance_matrix(data);

% Bus complex power injection
Si = - (data.bus(:, 3) + 1i * data.bus(:, 4));
genbuses = data.generator(:, 1);
Si(genbuses) = Si(genbuses) + (data.generator(:, 2) + 1i * data.generator(:, 3));
% Regulated voltages
Vr = zeros(data.nBuses, 1);
Vr(data.generator(:, 1)) = data.generator(:, 6);

% nInRow = 2 .* ones(2 * data.nBuses, 1);
nNonZeroInC = 0;
nNonZeroInD = 0;
nInc = zeros(data.nBuses, 1);
for i = 1:data.nBuses
    if data.bus(i, 2) == 1
        nInc(i) = nnz(data.powerSystemAC.nodalMatrix(i, :));
        nNonZeroInC = nNonZeroInC  + nInc(i) + 1;
%         nInRow(i) = incBranches + 1;
%         nInRow(i + data.nBuses) = incBranches + 1;
    elseif data.bus(i, 2) == 2
        nInc(i) = nnz(data.powerSystemAC.nodalMatrix(i, :));
        nNonZeroInD = nNonZeroInD + 2 * nInc(i) + 2;
%         nInRow(i) = 2 * incBranches;
    else
        nNonZeroInD = nNonZeroInD + 4;
    end
end

% memory allocation
rowsC = ones(nNonZeroInJ, 1);
colsC = zeros(nNonZeroInJ, 1);
valsC = zeros(nNonZeroInJ, 1);
rowsD = ones(nNonZeroInD, 1);
colsD = zeros(nNonZeroInD, 1);
valsD = zeros(nNonZeroInD, 1);
g = zeros(2 * data.nBuses, 1);

% initialization
x = [ data.bus(:, 8) .* exp(1i * data.bus(:, 9)), ...
      data.bus(:, 8) .* exp(1i * -data.bus(:, 9))  ];
iPQ = 0;
iPVTV = 0;
cntC = 1;
cntD = 1;


while iNo < solver.maxNumberOfIter
    for i = 1:data.nBuses
        if data.bus(i, 2) == 1
            iPQ = iPQ + 1;
            Ii = sum(data.powerSystemAC.nodalMatrix(i, :) .* x(1:nBuses));
            % calculate g
            g(i) = x(i) * conj(Ii)  - Si(i);
            g(i + data.nBuses) = conj(g(i)); 
            % calculate C
            nNew = nInc(i) + 1;
            idx = find(data.powerSystemAC.nodalMatrix(i, :));
            rowsC(cntC:cntC + nNew - 1) = iPQ;
            colsC(cntC:cntC + nNew - 1) = [i , data.nBuses + idx];
            valsC(cntC:cntC + nNew - 1) = [conj(Ii), x(i) .* data.powerSystemAC.nodalMatrix(i, :) ];
            cntC = cntC + nNew;
        elseif data.bus(i, 2) == 1
            iPVTV = iPVTV + 1;
            Ii = sum(data.powerSystemAC.nodalMatrix(i, :) .* x(1:nBuses));
            % calculate g
            g(i) = 0.5 * (x(i) * conj(Ii) + x(i + data.nBuses) * Ii) - real(Si(i));
            g(i + 1) = abs(x(i)) - Vr(i);
             % calculate D
            % active power
            nNew =  2 * nInc(i);
            idx = find(data.powerSystemAC.nodalMatrix(i, :));
            iIdx = find(idx == i);
            rowsD(cntD:cntD + nNew - 1) = iPVTV;
            valsS = [ zeros(1:iIdx - 1, 1), conj(Ii), zeros(iIdx + 1:idx - 1, 1),...
                      x(i) .* data.powerSystemAC.nodalMatrix(i, :) ];
            colsD(cntD:cntD + nNew - 1) = [idx , data.nBuses + idx];
            valsD(cntD:cntD + nNew - 1) = valsS + conj([valsS(nInc(i) + 1:2*nInc(i))...
                                          valsS(1:nInc(i))]);
            cntD = cntD + nNew;
            iPVTV = iPVTV + 1;
            % voltage magnitude
            rowsD(cntD:cntD + 1) = iPVTV;
            colsD(cntD:cntD + 1) = [i, data.nBuses + i];
            valsD(cntD:cntD + 1) = [x(data.nBuses + i), x(i)];
            cntD = cntD + 2;
        else 
            iPVTV = iPVTV + 1;
            % calculate g
            g(i) = abs(v(i)) - Vr(i);
            g(i) = imag(v(i));
            % calculate D
            rowsD(cntD:cntD + 3) = [ iPVTV, iPVTV, iPVTV + 1, iPVTV + 1 ];
            colsD(cntD:cntD + 1) = [ i, i + nBuses, i , i + nBuses ];
            valsD(cntD:cntD + 1) = [ 1/2, 1/2, -1i/2, 1i/2 ];
        end
    end
    rows = [ rowsC, rowsC + nPQ ];
    cols = [ colsC, mod((colsC + data.nBuses), 2 * data.nBuses) ];
    vals = [ valsC, conj(valsC) ];
end
end

