function [Xbnd,Tbnd,nodos_borde] = GetBoundary(X,intervalo,npx,npy,f1,f2,f3)
%devuelve 3 structs con las coordenadas, los elementos y los nodos del
%borde: por partes y completo
%hecho para el caso de una superficie "rectangular"
ax = intervalo(1);bx = intervalo(2);
ay = intervalo(3);by = intervalo(4);
xs = linspace(intervalo(1),intervalo(2),npx)'; 
ys = linspace(intervalo(3),intervalo(4),npy)';

%borde por partes
Xb1 = [f1(xs,ay*ones(size(xs))),f2(xs,ay*ones(size(xs))),f3(xs,ay*ones(size(xs)))];
[~,nodes_bnd1,~] = intersect(X,Xb1,'rows');
n_nodos_bnd1 = size(Xb1,1);
Tb1 = [(1:(n_nodos_bnd1-1))',(2:n_nodos_bnd1)']; 

Xb2 = [f1(bx*ones(size(ys)),ys),f2(bx*ones(size(ys)),ys),f3(bx*ones(size(ys)),ys)];
[~,nodes_bnd2,~] = intersect(X,Xb2,'rows');
n_nodos_bnd2 = size(Xb2,1);
Tb2 = [(1:(n_nodos_bnd2-1))',(2:n_nodos_bnd2)']; 

Xb3 = [f1(xs,by*ones(size(xs))),f2(xs,by*ones(size(xs))),f3(xs,by*ones(size(xs)))];
[~,nodes_bnd3,~] = intersect(X,Xb3,'rows');
n_nodos_bnd3 = size(Xb3,1);
Tb3 = [(1:(n_nodos_bnd3-1))',(2:n_nodos_bnd3)']; 

Xb4 = [f1(ax*ones(size(ys)),ys),f2(ax*ones(size(ys)),ys),f3(ax*ones(size(ys)),ys)];
[~,nodes_bnd4,~] = intersect(X,Xb4,'rows');
n_nodos_bnd4 = size(Xb4,1);
Tb4 = [(1:(n_nodos_bnd4-1))',(2:n_nodos_bnd4)']; 

%borde completo
Xb = [Xb1;Xb2;Xb3;Xb4]; n_nodos_bnd = size(Xb,1);
%orientamos el borde en el sentido contrario a las agujas del reloj
[~,nodes_bnd,~] = intersect(X,Xb,'rows');
Tb = [(1:(n_nodos_bnd-1))',(2:n_nodos_bnd)']; 
Tb = [Tb;[n_nodos_bnd,1]]; %cerramos el circulo

%construimos los structs
Xbnd=struct('Xb',Xb,'Xb1',Xb1,'Xb2',Xb2,'Xb3',Xb3,'Xb4',Xb4);
Tbnd=struct('Tb',Tb,'Tb1',Tb1,'Tb2',Tb2,'Tb3',Tb3,'Tb4',Tb4);
nodos_borde=struct('nodes_bnd',nodes_bnd,'nodes_bnd1',nodes_bnd1,...
    'nodes_bnd2',nodes_bnd2,'nodes_bnd3',nodes_bnd3,'nodes_bnd4',nodes_bnd4);