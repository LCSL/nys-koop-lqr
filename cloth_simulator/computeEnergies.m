function E = computeEnergies(phiPositions,phiVelocities,rho,Mlum,Fg,K)
n_frames = size(phiPositions,2);
E = zeros([1,n_frames]);
for i=1:n_frames
    phi = phiPositions{i}; 
    dphi = phiVelocities{i};
    E(i) = rho*0.5*dphi(:)'*Mlum*dphi(:) - phi(:)'*Fg + 0.5*phi(:)'*K*phi(:);
end 