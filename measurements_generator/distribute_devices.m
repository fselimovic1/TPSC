function [data] = distribute_devices(name, scada, pmu)
%---------------------------About Power System-----------------------------
mpc = ext2int(loadcase(name));
load(strcat('SG', name));

% SCADA
nScada = numel(scada.freqOfOccur);
scada.injP = [];
scada.injQ = [];
scada.flowP = [];
scada.flowQ = [];
scada.flowMagI = [];
scada.magV = [];

for i = 1:nScada
    if i == 1 && scada.freqOfOccur(1) ~= 0
        j = 1;
        if scada.freqOfOccur(1) < 0
            step = 1;
        else
            step = scada.freqOfOccur(1);
        end
        while j <= size(mpc.branch, 1)
            if scada.freqOfOccur(1) < 0 && mod(j, -scada.freqOfOccur(1)) == 0
                scada.flowP = [ scada.flowP;
                                mpc.branch(j, 2) i  -j scada.sd(1) scada.repRate ]; 
            end
            scada.flowP = [ scada.flowP;
                                mpc.branch(j, 1) i  j  scada.sd(1) scada.repRate ];
            j = j + step;
        end
    elseif i == 2 && scada.freqOfOccur(2) ~= 0
        j = 1;
        if scada.freqOfOccur(2) < 0
            step = 1;
        else
            step = scada.freqOfOccur(2);
        end
        while j <= size(mpc.branch, 1)
            if scada.freqOfOccur(2) < 0 && mod(j, -scada.freqOfOccur(2)) == 0
                scada.flowQ = [ scada.flowQ;
                                mpc.branch(j, 2) i  -j scada.sd(2) scada.repRate  ]; 
            end
            scada.flowQ = [ scada.flowQ;
                                mpc.branch(j, 1) i  j  scada.sd(2) scada.repRate  ];
            j = j + step;
        end
    elseif i == 3 && scada.freqOfOccur(3) ~= 0
        j = 1;
        while j <= size(mpc.bus, 1)
            scada.injP = [ scada.injP;
                           mpc.bus(j, 1) i  mpc.bus(j, 1) scada.sd(3) scada.repRate ];
            j = j + scada.freqOfOccur(3);
        end
    elseif i == 4 && scada.freqOfOccur(4) ~= 0
         j = 1;
        while j <= size(mpc.bus, 1)
            scada.injQ = [ scada.injQ;
                           mpc.bus(j, 1) i  mpc.bus(j, 1) scada.sd(4) scada.repRate ];
            j = j + scada.freqOfOccur(4);
        end
    elseif i == 5 && scada.freqOfOccur(5) ~= 0
        j = 1;
        if scada.freqOfOccur(5) < 0
            step = 1;
        else
            step = scada.freqOfOccur(5);
        end
        while j <= size(mpc.branch, 1)
            if scada.freqOfOccur(5) < 0 && mod(j, -scada.freqOfOccur(5)) == 0
                scada.flowMagI = [ scada.flowMagI;
                                    mpc.branch(j, 2) i  -j scada.sd(5) scada.repRate ]; 
            end
            scada.flowMagI = [ scada.flowMagI;
                               mpc.branch(j, 1) i  j scada.sd(5) scada.repRate ];
            j = j + step;
        end
    elseif i == 6 && scada.freqOfOccur(6) ~= 0
        j = 1;
        while j <= size(mpc.bus, 1)
            scada.magV = [ scada.magV;
                           mpc.bus(j, 1) i  mpc.bus(j, 1) scada.sd(6) scada.repRate];
            j = j + scada.freqOfOccur(6);
        end
    end
end

data.scada = [ scada.flowP; scada.flowQ; scada.injP; ...
    scada.injQ; scada.flowMagI; scada.magV ];

% Phasor Measurement Units
nPMUs = round(data.nBuses * pmu.dens / 100);
idx = randperm(data.nBuses);
idx = sort(idx(1:nPMUs)); 
randArray = make_rand_array(pmu.percDiffRepRates);
data.pmu = [ idx(1:nPMUs)', repmat([pmu.nCurrCh, pmu.sd], nPMUs, 1), ...
            randArray(randi(100, nPMUs, 1)) ];
% nPh = 3;
% pmu.ij = [];
% pmu.i = [];
% pmu.v = [];
% data.load = [];
% for i = 1:nPh
%     if i == 1 && template.pmu(1, 1) ~= 0
%         randArray = make_rand_array(template.pmu(1, 6:9));
%         j = 1;
%         if template.pmu(1, i) < 0
%             step = 1;
%         else
%             step = template.pmu(1, i);
%         end
%         while j <= size(mpc.branch, 1)
%             if template.pmu(1, i) < 0 && mod(j, -template.pmu(1, i)) == 0
%                 pmu.ij = [pmu.ij;
%                            1  -j template.pmu(1, 2:5) randArray(randi(100))]; 
%             end
%             pmu.ij = [ pmu.ij;
%                        1  j template.pmu(1, 2:5) randArray(randi(100))];
%             j = j + step;
%         end
%     elseif i == 2 && template.pmu(2, 1) ~= 0
%         randArray = make_rand_array(template.pmu(2, 6:9));
%         j = 1;
%         while j <= size(mpc.bus, 1)
%             pmu.ij = [pmu.ij;
%                        2  mpc.bus(j, 1) template.pmu(2, 2:5) randArray(randi(100))];
%             j = j + template.pmu(2, 1);
%         end
%     elseif i == 3 && template.pmu(3, 1) ~= 0 
%        randArray = make_rand_array(template.pmu(3, 6:9));
%         j = 1;
%         while j <= size(mpc.bus, 1)
%             pmu.ij = [pmu.ij;
%                        3  mpc.bus(j, 1) template.pmu(3, 2:5) randArray(randi(100))];
%             j = j + template.pmu(3, 1);
%         end
%     end   
% end
% data.pmu = [ pmu.ij; pmu.i; pmu.v ];



j = 1;
data.zi = [];
while j <= size(mpc.bus, 1)
     if mpc.bus(j, 2) == 1 && ...
             mpc.bus(j, 3) ==  0 && mpc.bus(j, 4) == 0
           if ~any(data.zi == j)
                data.zi = [data.zi; j];
           end
     end
    j = j + 1;
end
