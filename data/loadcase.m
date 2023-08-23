function data = loadcase(casename, varargin)
addpath("C:\Users\UIser\PowerSystemComputations\data\power_systems", ...
        "C:\Users\UIser\PowerSystemComputations\data\power_system_with_devices", ...
        "C:\Users\UIser\PowerSystemComputations\data\measurements");
if nargin == 1
    filename = strcat('TPSC', casename, '.mat');
elseif nargin == 2
    filename = strcat('TPSC', casename, 'D', '_', varargin{1}, '.mat');
elseif nargin == 3 
    if varargin{2} ~= 'M'
        ME = MException('MyComponent:notDefinedBehaviour', ...
                    'If entered, the thrid argument must be "M".');
        throw(ME); 
    else
        filename = strcat('TPSC', casename, 'M', '_', varargin{1}, '.mat');
    end
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
if nargin == 3
    data = measurements;
end
end

