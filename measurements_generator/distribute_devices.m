function [data] = distribute_devices(name, freqOfOccurSCADA,...
    freqOfOccurPMU, percDiffRepRates)
% SCADA - legacy measurements
% Measurement type (row number) - 
%               1 - Active power flow
%               2 - Reactive power flow
%               3 - Active power injection
%               4 - Rective power injection 
%               5 - Branch current magnitude
%               6 - Bus voltage magnitude   
%  1: - Frequency of occurrance 
%  2: - Standard deviation
%  3: - Reporting freqeuency
%
template.scada = [ 
                   freqOfOccurSCADA(1) 0.02 0.5
                   freqOfOccurSCADA(2) 0.02 0.5
                   freqOfOccurSCADA(3) 0.02 0.5
                   freqOfOccurSCADA(4) 0.02 0.5
                   freqOfOccurSCADA(5) 0.04 0.5
                   freqOfOccurSCADA(6) 0.002 0.5
                  ];
% MEASUREMENT DEVICES data format - SYNHROPHASOR MEASUREMENTS - WAMS
% 1: Measurement type - 
%               1 - Branch current phasor
%               2 - Injected current phasor
%               3 - Bus voltage phasor     
%
% 1: Step of appearance
% 2: Standard deviation (magnitude) [%]
% 3: Phase angle standard deviation [degrees]
% 4: Frequency standard deviationn [Hz] 
% 5: RoCoF standard deviation [Hz/s]
% 6: Reporting frequency 10 [Hz] - coefficient
% 7: Reporting frequency 25 [Hz] - coefficient
% 8: Reporting frequency 50 [Hz] - coefficient
% !!! 4: + 5: + 6: = 1 !!!
%
template.pmu = [
                  freqOfOccurPMU(1)  0.7  0.4  0.005  0  percDiffRepRates(1, :)
                  freqOfOccurPMU(2)  0.7  0.4  0.005  0  percDiffRepRates(2, :)
                  freqOfOccurPMU(3)  0.7  0.4  0.005  0  percDiffRepRates(3, :)
                 ];
%              
% % Load Dynamics
% % 1: Step of appearance
% % 2: Number of possible types
% template.load = [occursOfLoad 1];
  

%---------------------------About Power System-----------------------------
mpc = ext2int(loadcase(name));
load(strcat('SG', name));

% SCADA
nScada = size(template.scada, 1);
scada.injP = [];
scada.injQ = [];
scada.flowP = [];
scada.flowQ = [];
scada.flowMagI = [];
scada.magV = [];

for i = 1:nScada
    if i == 1 && template.scada(1, 1) ~= 0
        j = 1;
        if template.scada(1, 1) < 0
            step = 1;
        else
            step = template.scada(1, 1);
        end
        while j <= size(mpc.branch, 1)
            if template.scada(1, i) < 0 && mod(j, -template.scada(1, i)) == 0
                scada.flowP = [ scada.flowP;
                                1  -j template.scada(1, 2) template.scada(1, 3) ]; 
            end
            scada.flowP = [ scada.flowP;
                                1  j template.scada(1, 2) template.scada(1, 3) ];
            j = j + step;
        end
    elseif i == 2 && template.scada(2, 1) ~= 0
        j = 1;
        if template.scada(i, 1) < 0
            step = 1;
        else
            step = template.scada(i, 1);
        end
        while j <= size(mpc.branch, 1)
            if template.scada(2, i) < 0 && mod(j, -template.scada(2, i)) == 0
                scada.flowQ = [ scada.flowQ;
                                2  -j template.scada(1, 2) template.scada(1, 3) ]; 
            end
            scada.flowQ = [ scada.flowQ;
                                2  j template.scada(1, 2) template.scada(1, 3) ];
            j = j + step;
        end
    elseif i == 3 && template.scada(3, 1) ~= 0
        j = 1;
        while j <= size(mpc.bus, 1)
            scada.injP = [ scada.injP;
                           i  mpc.bus(j, 1) template.scada(i, 2)...
                           template.scada(i, 3)];
            j = j + template.scada(i, 1);
        end
    elseif i == 4 && template.scada(4, 1) ~= 0
         j = 1;
        while j <= size(mpc.bus, 1)
            scada.injQ = [ scada.injQ;
                           i  mpc.bus(j, 1) template.scada(i, 2)...
                           template.scada(i, 3)];
            j = j + template.scada(i, 1);
        end
    elseif i == 5 && template.scada(5, 1) ~= 0
        j = 1;
        if template.scada(i, 1) < 0
            step = 1;
        else
            step = template.scada(i, 1);
        end
        while j <= size(mpc.branch, 1)
            if template.scada(i, 1) < 0 && mod(j, - template.scada(i, 1)) == 0
                scada.flowMagI = [ scada.flowMagI;
                                    i  -j template.scada(i, 2) template.scada(i, 3) ]; 
            end
            scada.flowMagI = [ scada.flowMagI;
                               i  j template.scada(i, 2) template.scada(i, 3) ];
            j = j + step;
        end
    elseif i == 6 && template.scada(6, 1) ~= 0
        j = 1;
        while j <= size(mpc.bus, 1)
            scada.magV = [ scada.magV;
                           i  mpc.bus(j, 1) template.scada(i, 2) template.scada(i, 3)];
            j = j + template.scada(i, 1);
        end
    end
end

data.scada = [ scada.flowP; scada.flowQ; scada.injP; ...
    scada.injQ; scada.flowMagI; scada.magV ];

% Phasor Measurement Units;
nPh = size(template.pmu, 1);
pmu.ij = [];
pmu.i = [];
pmu.v = [];
data.zi = [];
data.load = [];
for i = 1:nPh
    if i == 1 && template.pmu(1, 1) ~= 0
        randArray = make_rand_array(template.pmu(1, 6:9));
        j = 1;
        if template.pmu(1, i) < 0
            step = 1;
        else
            step = template.pmu(1, i);
        end
        while j <= size(mpc.branch, 1)
            if template.pmu(1, i) < 0 && mod(j, -template.pmu(1, i)) == 0
                pmu.ij = [pmu.ij;
                           1  -j template.pmu(1, 2:5) randArray(randi(100))]; 
            end
            pmu.ij = [ pmu.ij;
                       1  j template.pmu(1, 2:5) randArray(randi(100))];
            j = j + step;
        end
    elseif i == 2 && template.pmu(2, 1) ~= 0
        randArray = make_rand_array(template.pmu(2, 6:9));
        j = 1;
        while j <= size(mpc.bus, 1)
            pmu.ij = [pmu.ij;
                       2  mpc.bus(j, 1) template.pmu(2, 2:5) randArray(randi(100))];
            j = j + template.pmu(2, 1);
        end
    elseif i == 3 && template.pmu(3, 1) ~= 0 
       randArray = make_rand_array(template.pmu(3, 6:9));
        j = 1;
        while j <= size(mpc.bus, 1)
            pmu.ij = [pmu.ij;
                       3  mpc.bus(j, 1) template.pmu(3, 2:5) randArray(randi(100))];
            j = j + template.pmu(3, 1);
        end
    end   
end
data.pmu = [pmu.ij; pmu.i; pmu.v];
j = 1;
while j <= size(mpc.bus, 1)
     if template.pmu(2, 1) == 0 && mpc.bus(j, 2) == 1 && ...
             mpc.bus(j, 3) ==  0 && mpc.bus(j, 4) == 0
           if ~any(data.zi == j)
                data.zi = [data.zi; j];
           end
     end
    j = j + 1;
end
