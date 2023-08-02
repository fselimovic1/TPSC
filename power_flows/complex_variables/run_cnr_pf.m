function results = run_cnr_pf(solver, data)
iNo = 1;
nPQ = sum(data.bus(:, 2) == 1);
nPV = sum(data.bus(:, 2) == 2);
nSlack = sum(data.bus(:, 2) == 3);

% initialization
x = [ data.bus(:, 8) .* exp(1i * data.bus(:, 9)); ...
      data.bus(:, 8) .* exp(1i * -data.bus(:, 9))  ];
  
% Compute admittance bus matrix
data.powerSystemAC = admittance_matrix(data);

% Bus complex power injection
Si = - (data.bus(:, 3) + 1i * data.bus(:, 4));
genbuses = data.generator(:, 1);
Si(genbuses) = Si(genbuses) + (data.generator(:, 2) + 1i * data.generator(:, 3));
% Regulated voltages
Vr = zeros(data.nBuses, 1);
Vr(data.generator(:, 1)) = data.generator(:, 6);

nNonZeroInC = 0;
nNonZeroInD = 0;
slack = [];
nInc = zeros(data.nBuses, 1);
for i = 1:data.nBuses
    if data.bus(i, 2) == 1
        nInc(i) = nnz(data.powerSystemAC.nodalMatrix(i, :));
        nNonZeroInC = nNonZeroInC  + nInc(i) + 1;
    elseif data.bus(i, 2) == 2
        nInc(i) = nnz(data.powerSystemAC.nodalMatrix(i, :));
        nNonZeroInD = nNonZeroInD + 2 * nInc(i) + 2;
    else
        nNonZeroInD = nNonZeroInD + 4;
    end
end

% memory allocation
rowsC = zeros(nNonZeroInC, 1);
colsC = zeros(nNonZeroInC, 1);
rowsD = zeros(nNonZeroInD, 1);
colsD = zeros(nNonZeroInD, 1);
valsC = zeros(nNonZeroInC, 1);
valsD = zeros(nNonZeroInD, 1);
g = zeros(2 * data.nBuses, 1);

iPQ = 0;
iPVTV = 0;
cntC = 1;
cntD = 1;
for i = 1:data.nBuses
    if data.bus(i, 2) == 1
        iPQ = iPQ + 1;
        Ii = sum(data.powerSystemAC.nodalMatrix(i, :).' .* x(1:data. nBuses));
        % calculate g
        g(i) = x(i) * conj(Ii)  - Si(i);
        g(i + data.nBuses) = conj(g(i)); 
        % calculate C
        nNew = nInc(i) + 1;
        idx = find(data.powerSystemAC.nodalMatrix(i, :));
        rowsC(cntC:cntC + nNew - 1) = iPQ;
        colsC(cntC:cntC + nNew - 1) = [i , data.nBuses + idx];
        valsC(cntC:cntC + nNew - 1) = [ conj(Ii), nonzeros(x(i) .* data.powerSystemAC.nodalMatrix(i, :))' ];
        cntC = cntC + nNew;
    elseif data.bus(i, 2) == 2
        iPVTV = iPVTV + 1;
        Ii = sum(data.powerSystemAC.nodalMatrix(i, :).' .* x(1:data.nBuses));
        % calculate g
        g(i) = 0.5 * (x(i) * conj(Ii) + x(i + data.nBuses) * Ii) - real(Si(i));
        g(i + 1) = abs(x(i)) - Vr(i);
        % calculate D
        % active power
        nNew =  2 * nInc(i);
        idx = find(data.powerSystemAC.nodalMatrix(i, :));
        iIdx = find(idx == i);
        rowsD(cntD:cntD + nNew - 1) = iPVTV;
        valsS = [ zeros(1, iIdx - 1), conj(Ii), zeros(1, nInc(i) - iIdx),...
                  nonzeros(x(i) .* data.powerSystemAC.nodalMatrix(i, :)).' ];
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
         % calculate g
         g(i) = abs(x(i)) - Vr(i);
         g(i) = imag(x(i));
         % calculate D
         rowsD(cntD:cntD + 3) = [ iPVTV + 1, iPVTV + 1, iPVTV + 2, iPVTV + 2 ];
         colsD(cntD:cntD + 3) = [ i, i + data.nBuses, i , i + data.nBuses ];
         valsD(cntD:cntD + 3) = [ 1/2, 1/2, -1i/2, 1i/2 ];
         iPVTV = iPVTV + 2;
         cntD = cntD + 4;
    end
end
rows = [ rowsC; rowsC + nPQ; rowsD + 2 * nPQ ];
cols = [ colsC; mod((colsC + data.nBuses), 2 * data.nBuses); colsD ];
cols(cols == 0) = 2 * data.nBuses;
vals = [ valsC; conj(valsC); valsD ];

J = sparse(rows, cols, vals, 2 * data.nBuses, 2 * data.nBuses);

while iNo < solver.maxNumberOfIter
    dx = - J \ g;
    x = x + dx;
    
    cntC = 1;
    cntD = 1;
    for i = 1:data.nBuses
        if data.bus(i, 2) == 1
            Ii = sum(data.powerSystemAC.nodalMatrix(i, :)' .* x(1:nBuses));
            % calculate g
            g(i) = x(i) * conj(Ii)  - Si(i);
            g(i + data.nBuses) = conj(g(i)); 
            % calculate C
            nNew = nInc(i) + 1;
            valsC(cntC:cntC + nNew - 1) = [ conj(Ii), x(i) .* data.powerSystemAC.nodalMatrix(i, :) ];
            cntC = cntC + nNew;
        elseif data.bus(i, 2) == 1
            Ii = sum(data.powerSystemAC.nodalMatrix(i, :)' .* x(1:nBuses));
            % calculate g
            g(i) = 0.5 * (x(i) * conj(Ii) + x(i + data.nBuses) * Ii) - real(Si(i));
            g(i + 1) = abs(x(i)) - Vr(i);
            % calculate D
            % active power
            nNew =  2 * nInc(i);
            idx = find(data.powerSystemAC.nodalMatrix(i, :));
            iIdx = find(idx == i);
            valsS = [ zeros(1:iIdx - 1, 1), conj(Ii), zeros(iIdx + 1:idx - 1, 1),...
                      x(i) .* data.powerSystemAC.nodalMatrix(i, :) ];
            valsD(cntD:cntD + nNew - 1) = valsS + conj([valsS(nInc(i) + 1:2*nInc(i))...
                                          valsS(1:nInc(i))]);
            cntD = cntD + nNew;
            % voltage magnitude
            valsD(cntD:cntD + 1) = [x(data.nBuses + i), x(i)];
            cntD = cntD + 2;
        else 
            % calculate g
            g(i) = abs(v(i)) - Vr(i);
            g(i) = imag(v(i));
            % calculate D
            valsD(cntD:cntD + 1) = [ 1/2, 1/2, -1i/2, 1i/2 ];
        end
    end
    
    
end
end

