function [results] = run_power_flows(pfsettings, data, varargin)
if nargin == 3
    pfsettings.postprocess = 1;
    if strcmp('complex', pfsettings.domain)
        results = run_cnr_pf(pfsettings, data, varargin{1});
    end
else
    if strcmp('complex', pfsettings.domain)
        results = run_cnr_pf(pfsettings, data);
    end
end
end
