clear all
close all
clc

load('LOG/GENERATED_DATA/TruePhiPositions.mat');
load('LOG/GENERATED_DATA/PredPhiPositions.mat');

flag = 0;
if flag == 0
    PhiPositions = TruePhiPositions;
else
    PhiPositions = PredPhiPositions;
end

ax = -0.5; bx = 0.5; ay = -0.5; by = 0.5; 
f1 = @(x,y) x; f2 = @(x,y) 0.75*sqrt(1 - x.^2) +0.1; f3 = @(x,y) y;
npx = 8; npy = 8; n_nodos = npx*npy; %numero de puntos
[X,T] = CreateMesh([ax,bx,ay,by],npx,npy,f1,f2,f3); 
[Xb,Tb,nodos_borde] = GetBoundary(X,[ax,bx,ay,by],npx,npy,f1,f2,f3); 

peli = 'video'; vel_vid = 1; 
VID = makeVideo(PhiPositions,T,nodos_borde,peli,vel_vid,0.05)