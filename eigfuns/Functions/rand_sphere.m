% Generates uniform samples over the unit ball

function x_samp = rand_sphere(n,n_samp)

if(~exist('n_samp','var'))
    n_samp = 1;
end

x_samp = [];
while(size(x_samp,2) < n_samp)
    x_samp = [x_samp, 2*rand(n,n_samp)-1];
    x_samp = x_samp(:, sum(x_samp.^2) <= 1 );
end
x_samp = x_samp(:,1:n_samp);

end

