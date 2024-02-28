function [Cphi,DCphi] = fun_C(phi,C,Mcons,A_b,n_nodos,n_conds)
%devuelve la funcion de restricciones y su gradiente 
%evaluada en la superficie actual
conds = reshape(C',[n_conds*n_nodos,n_nodos])*reshape(phi,[n_nodos,3]);
%Jacobiano
conds_pos = reshape(conds(:,1:3),[n_conds,3*n_nodos]); 
grad = conds_pos;
%Valor de la funcion
Cphi = sparse([grad;A_b])*reshape(phi,[3*n_nodos,1]);
DCphi = sparse([2*grad;A_b]);
