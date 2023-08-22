function data = loadcase(casename, varargin)
addpath("C:\Users\UIser\PowerSystemComputations\data\power_systems");
addpath("C:\Users\UIser\PowerSystemComputations\data\power_system_with_devices");
if nargin == 1
    filename = strcat('TPSC', casename, '.mat');
elseif nargin == 2
    filename = strcat('TPSC', casename, 'D', '_', varargin{1}, '.mat');
else
    ME = MException('MyComponent:notDefinedBehaviour', ...
                    'Inappropriate number of input arguments for the function "loadcase".');
    throw(ME); 
end
if ~exist(filename, 'file')
	ME = MException('MyComponent:notDefinedBehaviour', ...
                     'Specifed file does not exist.');
	throw(ME); 
else
	load(filename);
end
end

