function [ measurements ] = generatemeasurements(data)
% number of devices
nPMU = size(data.pmu, 1);
nSCADA  = size(data.scada, 1);

% In tracking mode, calculating frequency is equal to the maximum reporting     
% frequency of devices embedded in the grid.
if strcmp(data.mode, 'tracking')
    if ~nPMU
        calcfreq = max(data.pmu(:, 7)); 
    else
        calcfreq = max(data.scada(:, 5));
    end
    tstamps = data.t * calcfreq;
    if round(tstamps) ~= tstamps
        tstamps = round(tstamps); 
    else
        tstamps = tstamps + 1;
    end
    % number of measurements in the interval
    mPMU = 0;
    mSCADA = 0;
    for i = 1:nPMU
        mPMU = mPMU + fix(tstamps * data.pmu(i, 7) * calcfreq);
    end
    for i = 1:nPMU
        mSCADA = mSCADA + fix(tstamps * data.scada(i, 5) * calcfreq);
    end
end
end

