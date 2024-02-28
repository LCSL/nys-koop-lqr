clear all
close all
clc

n_points = 8;

true_data = csvread('LOG/REAL_DATA/test_samples_14.csv');
test_data = csvread('LOG/REAL_DATA/test_samples_hat_14.csv');

[n_step, ~] = size(true_data);

TruePhiPositions = {};
for t=1:n_step
    mesh = zeros(n_points^2,3);
    mesh(:,1) = true_data(t,1:n_points^2); % X
    mesh(:,2) = true_data(t,n_points^2+1:2*n_points^2); % Y
    mesh(:,3) = true_data(t,2*n_points^2+1:end); % Z
    TruePhiPositions{t} = mesh;
end

PredPhiPositions = {};
for t=1:n_step
    mesh = zeros(n_points^2,3);
    mesh(:,1) = test_data(t,1:n_points^2); % X
    mesh(:,2) = test_data(t,n_points^2+1:2*n_points^2); % Y
    mesh(:,3) = test_data(t,2*n_points^2+1:end); % Z
    PredPhiPositions{t} = mesh;
end

flag = 0;
if flag == 0
    PhiPositions = TruePhiPositions;
else
    PhiPositions = PredPhiPositions;
end


ax = -0.25; bx = 0.5; ay = -0.5; by = 0.5; 
f1 = @(x,y) x; f2 = @(x,y) 0.75*sqrt(1 - x.^2) +0.1; f3 = @(x,y) y;
npx = 8; npy = 8; n_nodos = npx*npy; %numero de puntos
[X,T] = CreateMesh([ax,bx,ay,by],npx,npy,f1,f2,f3); 
[Xb,Tb,nodos_borde] = GetBoundary(X,[ax,bx,ay,by],npx,npy,f1,f2,f3); 

peli = 'video'; vel_vid = 1; 
VID = makeVideo(PhiPositions,T,nodos_borde,peli,vel_vid,0.05)