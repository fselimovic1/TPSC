function freq_profile = get_predefined_curve(fM, tSE, tS, fn, fType, fSubType, fMin)
if fType == 1
    % to be done
end
if fType == 2
    t = 0:1/fM:tSE - 1/fM;
    freq_profile = fn + (fn - fMin) * sin(2 * pi * t); 
end
if fType == 3
    if fSubType == 1
        tDyn = 0.2 * tSE;
        freq_profile = curve1(fM, fn, fMin, tS, tSE - tS, tDyn);
    elseif fSubType == 2
        tE = 0;
        freq_profile = curve2(fM, tS, tE, tSE - tS - tE, fn, fMin);
    end
%     fc = 450; % cutoff frequency
%     fs = 1000; % sample rate
%     % Filter to achieve the smootheness
%     [b, a] = butter(6, fc/(fs/2), 'low'); % 6th order low-pass Butterworth filter
%     % Apply the filter
%     freq_profile = filter(b,a,freq_profile);
%     freq_profile(1:0.8 * fM) = fn;
end
end


