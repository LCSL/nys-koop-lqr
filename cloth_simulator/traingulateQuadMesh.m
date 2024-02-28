function [Xtri,Ttri] =  traingulateQuadMesh(X,T)
n_elem = size(T,1);
k1 = T(:,1); k2 = T(:,2); k3 = T(:,3); k4 = T(:,4); k5 = size(X,1)+[1:n_elem]';
Xint= (X(k1,:) + X(k2,:) + X(k3,:) + X(k4,:))/4;
Xtri = [X;Xint];
Ttri = [k1 k2 k5; k2 k3 k5; k3 k4 k5; k4 k1 k5];
end