function [super,curve] = updatePlot(super,curve,X,T,nodos_borde)
phi = X;
%gamma = X(nodos_borde.nodes_bnd,:);
k1 = T(:,1); k2 = T(:,2); k3 = T(:,3); k4 = T(:,4);
phi_int= (phi(k1,:) + phi(k2,:) + phi(k3,:) + phi(k4,:))/4;
super.Vertices = [phi;phi_int];
%curve.XData = gamma(:,1);curve.YData = gamma(:,2);curve.ZData = gamma(:,3);
end