function [X,T] = CreateMesh(intervalo,npx,npy,f1,f2,f3)
% Crea un mallado y una cuadrangulacion dado un rectangulo [ax,bx]x[ay,by]
% y una funcion altura u, e.g. u(x,y) = sqrt(1-x^2)
% dim X = n_puntos x 3 (coordenadas en el espacio)
% dim T = n_elementos x 4 (indices de los nodos)

% Allocate space for the nodal coordinates matrix
X = zeros((npx)*(npy),3);
xs = linspace(intervalo(1),intervalo(2),npx)'; 
unos = ones(npx,1);
% Nodes' coordinates
yys = linspace(intervalo(3),intervalo(4),npy);
for i=1:npy
    ys = yys(i)*unos; 
    posi = (i-1)*(npx)+1:i*(npx); 
    X(posi,:)=[f1(xs,ys),f2(xs,ys),f3(xs,ys)];
end

%elementos (cuadrilateros)
nx = npx-1; ny = npy-1;
T = zeros(nx*ny,4);
for a=1:ny
    for b=1:nx
        ielem = (a-1)*nx+b;
        inode = (a-1)*(npx)+b;
        T(ielem,:) = [inode   inode+1   inode+npx+1   inode+npx];
    end   
end

