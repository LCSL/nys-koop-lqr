function error = computeAreaError(phiPositions,T,theReferenceElement)
A0 = computeArea(phiPositions{1},T,theReferenceElement);
n_frames = size(phiPositions,2);
error = zeros([1,n_frames]);
for i=2:n_frames
    Af = computeArea(phiPositions{i},T,theReferenceElement);
    error(i) = 100*abs(sum(A0)-sum(Af))/sum(A0);
end    