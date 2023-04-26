function [ data ] = compute_branch_current_coeffs(data)
data.powerSystemAC.A = complex(zeros(data.nBranches, 1));
data.powerSystemAC.B = complex(zeros(data.nBranches, 1));
data.powerSystemAC.C = complex(zeros(data.nBranches, 1));
data.powerSystemAC.D = complex(zeros(data.nBranches, 1));
for i = 1:data.nBranches
    % if the branch is IN-SERVICE:
    if data.branch(i, 11)
        yij = 1 / (data.branch(i, 3) + 1i * data.branch(i, 4));
        gij = real(yij); bij = imag(yij);
        bs = data.branch(i, 4) / 2;
        data.powerSystemAC.A(i) = gij^2 + (bij + bs)^2;
        data.powerSystemAC.B(i) = gij^2 + bij^2;
        data.powerSystemAC.C(i) = gij^2 + bij(bij + bs);
        data.powerSystemAC.D(i) = gij * bs;
    end
end
end

