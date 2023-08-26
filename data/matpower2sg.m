% Builds the power system data form Matpower7.1. and process data to make
% it usable.
clc
clear

% insert the name of Matpower case 
caseName = 'case13659pegase';

data = eval(caseName);

% DATA BUS
data.bus(:, 3:6) = data.bus(:, 3:6)./data.baseMVA;
data.bus(:, 9) = data.bus(:, 9) .* pi/180;
data.bus(:, 14:15) = zeros(size(data.bus, 1), 2);

% DATA GENERATOR
data.generator = data.gen(:, 1:10);
data.generator(:, [2:5, 9, 10]) = data.generator(:, [2:5, 9, 10])./data.baseMVA;
data = rmfield(data, 'gen');

% for i = 1:size(data.generator, 1)
%     data.bus(data.generator(i, 1), 14) = data.bus(data.generator(i, 1), 14) + data.generator(i, 2);
%     data.bus(data.generator(i, 1), 15) = data.bus(data.generator(i, 1), 15) + data.generator(i, 3);
% end
% DATA BRANCH
data.branch(:, 6:8) = data.branch(:, 6:8)./data.baseMVA;
data.branch(:, [10, 12, 13]) = data.branch(:, [10, 12, 13]) .* pi/180;

% NEW
if sum(data.bus(:, 2) == 3) > 1
    assert(false)
end
data.fn = 50;
data.case = caseName;

home = getenv('USERPROFILE');
path = strcat(home, '\PowerSystemComputations\data\power_systems\TPSC', caseName);
save(path,'data')

