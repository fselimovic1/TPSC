function [h_i, H_i] = c_slack_angle(nBuses, nSlack, state)
% measurement function
h_i = 1i/2 * (state(nSlack + nBuses) - state(nSlack));

% measurements jacobian matrix
H_i = sparse([1, 1], [nSlack, nBuses + nSlack], [-1i/2, 1i/2], 1, 2 * nBuses);
