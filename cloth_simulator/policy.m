function out = policy(state,err,lengthscales,centers,weights,u_max)
%POLICY Sum of gaussians policy
%   Detailed explanation goes here

num_basis = length(weights);

% extended state
ext_state = [state; err];

% normalize states and centers
norm_state = ext_state./lengthscales;
norm_centers = (centers./repmat(lengthscales,1,num_basis))';

% get the square distance
dist = sum(norm_state.^2);
dist =  dist + sum(norm_centers.^2,2);
dist = dist - 2*norm_centers*norm_state;

% weighted sum of gaussians
delta_u = weights'*exp(-dist);

% squash
delta_u = u_max*tanh(delta_u/u_max);

out = [delta_u; delta_u];
end

