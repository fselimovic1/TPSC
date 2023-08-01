function results = run_cnr_pf(solver, data)
% Compute admittance bus matrix
data.powerSystemAC = admittance_matrix(data);

% Bus complex power injection
Sinj = - (data.bus(:, 3) + 1i * data.bus(:, 4));
genbuses = data.generator(:, 1);
Sinj(genbuses) = Sinj(genbuses) + (data.generator(:, 2) + 1i * data.generator(:, 3));

end

