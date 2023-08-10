function [results] = run_power_flows(pfsettings, data, varargin)
if strcmp('complex', pfsettings.domain)
    results = run_cnr_pf(pfsettings, data, varargin);
end
end
