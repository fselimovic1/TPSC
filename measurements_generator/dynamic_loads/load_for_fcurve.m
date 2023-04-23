function loadC = load_for_fcurve(fSubType, fM, fn, f)
loadC = 1 - f * 0;
loadC(fM + 1:end) = 1.1;
end

