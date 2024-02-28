function L = long_aristas(puntos,quad)
n = size(quad,1);
L = 0;
for i = 1:n 
    t = quad(i,:);
    a = puntos(t(1),:); b = puntos(t(2),:); 
    L = L + norm(b-a);
end


