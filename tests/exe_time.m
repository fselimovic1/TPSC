clc
clear 

%-------------------------- Average execution time ------------------------
%--------------------------------------------------------------------------

nTrails = 50;
times = zeros(nTrails, 1);
%------------------------- Power System Case ------------------------------
casename = 'case1888rte';
vrs = 'A';
%--------------------------------------------------------------------------

%------------------- State estimaion solver - Settings --------------------
sesettings.domain = 'complex';
sesettings.method = 'ecgn_sse';
sesettings.mweights = [ "pmuscadaratio", 5 ];
sesettings.virtual = 0;
sesettings.flatStart = 0;
sesettings.maxNumberOfIter = 50;
sesettings.eps = 1e-6;
sesettings.showresults = 0;
%--------------------------------------------------------------------------

%------------------------ Run state estimation ----------------------------
for i = 1:nTrails
    results = runsse(casename, vrs, sesettings);
    times(i) = results.algtime;
end
clc
fprintf("Average execution time is: %.2f.\n", mean(times) * 1000);
%--------------------------------------------------------------------------







