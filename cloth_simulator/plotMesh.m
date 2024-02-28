function [super,curve] = plotMesh(X,T,nodos_borde)
figure(1)
super = trisurf(T,X(:,1),X(:,2),X(:,3),'edgecolor', 'k','facecolor','w');
% super = trisurf(T,X(:,1),X(:,2),X(:,3),'edgecolor', 'k','facecolor',[0.4471    0.9020    0.3451]);

hold on
%gamma0 = X(nodos_borde.nodes_bnd,:);
curve = [];% plot3(gamma0(:,1),gamma0(:,2),gamma0(:,3),'.','color','k','MarkerSize',10);
% light               % add a light
% lighting gouraud    % preferred lighting for a curved surface
%ylim([-0.65,0.0]); xlim([-0.35,0.3]); zlim([-0.15,0.5]) % REAL CLOTH LIM
ylim([-2,2]); xlim([-1.5,1.5]); zlim([-1.5,1.5]) % SIMULATION LIM
camzoom(1.5)   % zoom into scene
axis off
% axis tight
% axis equal
% set axis equal and remove axis
xlabel('X'); ylabel('Y');zlabel('Z');