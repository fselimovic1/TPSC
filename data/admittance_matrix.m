function powerSystemAC = admittance_matrix(data, varargin)
complexVector = complex(zeros(data.nBranches, 1));

powerSystemAC.transformerRatio = complexVector;
powerSystemAC.admittance = complexVector;
powerSystemAC.nodalToTo = complexVector;
powerSystemAC.nodalFromFrom = complexVector;
powerSystemAC.nodalFromTo = complexVector;
powerSystemAC.nodalToFrom = complexVector;
powerSystemAC.from = data.branch(:, 1);
powerSystemAC.to = data.branch(:, 2);
nodalDiagonals = complex(zeros(data.nBuses, 1));
if nargin == 2
     for i =  1:data.nBranches
        % if the branch is IN-SERVICE:
        if data.branch(i, 11)
            powerSystemAC.admittance(i) = 1 / (data.branch(i, 3) + 1i * data.branch(i, 4) * varargin{1}.f / data.fn);
        
            if  ~data.branch(i, 9)
                powerSystemAC.transformerRatio(i) = exp(1i * data.branch(i, 10));
            else
                powerSystemAC.transformerRatio(i) = data.branch(i, 9) * exp(1i * data.branch(i, 10));
            end
            transformerRatioConj = conj(powerSystemAC.transformerRatio(i));
        
            % Branch PI model elements:
            powerSystemAC.nodalToTo(i) =  powerSystemAC.admittance(i) + 1i / 2 * (data.branch(i, 5) * (data.fn / varargin{1}.f));
            powerSystemAC.nodalFromFrom(i) =  powerSystemAC.nodalToTo(i) / (transformerRatioConj * powerSystemAC.transformerRatio(i));
            powerSystemAC.nodalFromTo(i) = -powerSystemAC.admittance(i) / transformerRatioConj;
            powerSystemAC.nodalToFrom(i) = -powerSystemAC.admittance(i) / powerSystemAC.transformerRatio(i);
        
            nodalDiagonals(data.branch(i, 1)) =  nodalDiagonals(data.branch(i, 1)) + powerSystemAC.nodalFromFrom(i);
            nodalDiagonals(data.branch(i, 2)) =  nodalDiagonals(data.branch(i, 2)) + powerSystemAC.nodalToTo(i);
        end
    end

    for i = 1:data.nBuses
        nodalDiagonals(i) = nodalDiagonals(i) + data.bus(i, 5) + 1i * data.bus(i, 6) * (data.fn / varargin{1}.f) ;
    end
else
    for i = 1:data.nBranches
        % if the branch is IN-SERVICE:
        if data.branch(i, 11)
            powerSystemAC.admittance(i) = 1 / (data.branch(i, 3) + 1i * data.branch(i, 4));
        
            if  ~data.branch(i, 9)
                powerSystemAC.transformerRatio(i) = exp(1i * data.branch(i, 10));
            else
                powerSystemAC.transformerRatio(i) = data.branch(i, 9) * exp(1i * data.branch(i, 10));
            end
            transformerRatioConj = conj(powerSystemAC.transformerRatio(i));
        
            % Branch PI model elements:
            powerSystemAC.nodalToTo(i) =  powerSystemAC.admittance(i) + 1i / 2 * data.branch(i, 5);
            powerSystemAC.nodalFromFrom(i) =  powerSystemAC.nodalToTo(i) / (transformerRatioConj * powerSystemAC.transformerRatio(i));
            powerSystemAC.nodalFromTo(i) = -powerSystemAC.admittance(i) / transformerRatioConj;
            powerSystemAC.nodalToFrom(i) = -powerSystemAC.admittance(i) / powerSystemAC.transformerRatio(i);
        
            nodalDiagonals(data.branch(i, 1)) =  nodalDiagonals(data.branch(i, 1)) + powerSystemAC.nodalFromFrom(i);
            nodalDiagonals(data.branch(i, 2)) =  nodalDiagonals(data.branch(i, 2)) + powerSystemAC.nodalToTo(i);
        end
    end

    for i = 1:data.nBuses
        nodalDiagonals(i) = nodalDiagonals(i) + data.bus(i, 5) + 1i * data.bus(i, 6);
    end
end
% Compose sparse matrices:
powerSystemAC.nodalMatrix = sparse([(1:data.nBuses)'; powerSystemAC.from;...
                                    powerSystemAC.to], [(1:data.nBuses)'; ...
                                    powerSystemAC.to; powerSystemAC.from], ...
                            [ nodalDiagonals; powerSystemAC.nodalFromTo; powerSystemAC.nodalToFrom], ...
                             data.nBuses, data.nBuses);
powerSystemAC.nodalMatrixTranspose = transpose(powerSystemAC.nodalMatrix);  
powerSystemAC.nodalMatrixOffDiag = sparse([powerSystemAC.from; powerSystemAC.to], ...
                                [ powerSystemAC.to; powerSystemAC.from], ...
                                [ powerSystemAC.nodalFromTo; powerSystemAC.nodalToFrom], ...
                                data.nBuses, data.nBuses);
end


