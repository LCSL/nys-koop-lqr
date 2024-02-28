clear all
close all
clc
load('LOG/log_gpdm.mat')

Y = csvread('LOG/Y_rollout.csv')';
U = csvread('LOG/U_rollout.csv')';
[state_dim, ~] = size(Y);
[input_dim, num_steps] = size(U); 
num_free_nodes = state_dim/3;
num_input_nodes = input_dim/3;
num_nodes = num_free_nodes+num_input_nodes;

free_nodes = linspace(1,num_nodes,num_nodes);
free_nodes(controlados) = [];

for t=1:num_steps
    input_nodes = reshape(U(:,t), [3,num_input_nodes])';
    state_nodes = reshape(Y(:,t), [3,num_free_nodes])';
    newPhiPositions{t} = zeros(num_nodes,3);
    newPhiPositions{t}(controlados,:) = input_nodes;
    newPhiPositions{t}(free_nodes,:) = state_nodes;
    
end
save('LOG/ControlPhiPositions.mat', 'newPhiPositions');



