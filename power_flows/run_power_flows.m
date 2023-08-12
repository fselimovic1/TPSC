function [results] = run_power_flows(pfsettings, data)
if strcmp('complex', pfsettings.domain)
    results = run_cnr_pf(pfsettings, data);
end
end
