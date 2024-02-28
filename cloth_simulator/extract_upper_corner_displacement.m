clear all
close all
clc
path = "/Users/edoardocaldarelli/Documents/phd/MPC/SIMULADORforCGPDM/LOG/CONTROL_DATA/";
npx = 8; npy = 8; n_nodos = npx*npy;
controlados = [(n_nodos - npx)* 3 + 1, (n_nodos - npx)* 3 + 2, (n_nodos - npx)* 3 + 3, n_nodos * 3 - 2, n_nodos * 3 - 1, n_nodos * 3];
for i = (0:29)
    curr_state = readmatrix(path + "state_samples_cloth_swing_" + i + ".csv");
    delta_states = diff(curr_state);
    delta_u = delta_states(:, controlados);
end
