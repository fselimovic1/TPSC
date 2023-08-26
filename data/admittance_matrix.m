function powerSystemAC = admittance_matrix(num, bus, branchi, branchj, resistance, reactance, charging, transturns, transphase, branchstatus, Gshunt, Bshunt)
% allocate memory
complexVector = complex(zeros(num.branch, 1));
powerSystemAC.transformerRatio = complexVector;
powerSystemAC.admittance = complexVector;
powerSystemAC.nodalToTo = complexVector;
powerSystemAC.nodalFromFrom = complexVector;
powerSystemAC.nodalFromTo = complexVector;
powerSystemAC.nodalToFrom = complexVector;
nodalDiagonals = complex(zeros(num.bus, 1));

powerSystemAC.admittance(branchstatus) = 1 ./ (resistance(branchstatus) + 1i .* reactance(branchstatus));
isLine = ~transturns;
powerSystemAC.transformerRatio(isLine) = exp(1i .* transphase(isLine));
powerSystemAC.transformerRatio(~isLine) = transturns(~isLine) .* exp(1i * transphase(~isLine));
transformerRatioConj = conj(powerSystemAC.transformerRatio);
% Branch PI model elements:
powerSystemAC.nodalToTo =  powerSystemAC.admittance + 1i / 2 .* charging;
powerSystemAC.nodalFromFrom =  powerSystemAC.nodalToTo ./ (transformerRatioConj .* powerSystemAC.transformerRatio);
powerSystemAC.nodalFromTo = -powerSystemAC.admittance ./ transformerRatioConj;
powerSystemAC.nodalToFrom = -powerSystemAC.admittance ./ powerSystemAC.transformerRatio;
nodalDiagonals = nodalDiagonals + accumarray([branchi; branchj ], [powerSystemAC.nodalFromFrom; powerSystemAC.nodalToTo ], [num.bus, 1] );

% + Shunts
nodalDiagonals = nodalDiagonals + Gshunt + 1i .* Bshunt;

% Compose sparse matrices:
powerSystemAC.nodalMatrix = sparse([(bus)'; branchi;...
                                    branchj], [(bus)'; ...
                                   branchj; branchi], ...
                            [ nodalDiagonals; powerSystemAC.nodalFromTo; powerSystemAC.nodalToFrom], ...
                             num.bus, num.bus);
powerSystemAC.nodalMatrixTranspose = transpose(powerSystemAC.nodalMatrix);  
end



