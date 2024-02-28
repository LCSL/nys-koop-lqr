function error = computeEdgesError(phiPositions,Tb)
L0 = long_aristas(phiPositions{1},Tb);
n_frames = size(phiPositions,2);
error = zeros([1,n_frames]);
for i=2:n_frames
    Lf = long_aristas(phiPositions{i},Tb);
    error(i) = 100*abs(L0-Lf)./L0;
end 