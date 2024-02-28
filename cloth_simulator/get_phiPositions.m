clear all
close all
clc

num = 9;
fname = strcat('LOG/GENERATED_DATA/test_samples_',num2str(num),'.csv');
Ytest = csvread(fname)';

[num_states, num_steps] = size(Ytest); 
num_nodes = num_states/3;
% Build State samples
for t=1:num_steps
    TruePhiPositions{t} = reshape(Ytest(:,t), [3,num_nodes])';
end

save('LOG/GENERATED_DATA/TruePhiPositions.mat', 'TruePhiPositions');

fname = strcat('LOG/GENERATED_DATA/test_samples_hat_',num2str(num),'.csv');
Yhat = csvread(fname)';

[num_states, num_steps] = size(Yhat); 
num_nodes = num_states/3;
% Build State samples
for t=1:num_steps
    PredPhiPositions{t} = reshape(Yhat(:,t), [3,num_nodes])';
end

save('LOG/GENERATED_DATA/PredPhiPositions.mat', 'PredPhiPositions');