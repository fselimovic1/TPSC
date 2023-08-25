function [ results, data ] = runpf(casename, pfsettings)
data = loadcase(casename);
if strcmp('complex', pfsettings.domain)
	results = run_cnr_pf(pfsettings, data);
    results_pf(data, pfsettings, results);
end
end
