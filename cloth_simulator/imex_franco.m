function [phi,dphi] = imex_franco(dt,phi0,dphi0,phi1,dphi1,...
                                   Fg,rho,A_b,Minv,Mcons,K,D,C,...
                                   n_conds,n_nodos,C0,u)
%hacemos un reshape
phi0 = reshape(phi0,[3*n_nodos,1]);
phi1 = reshape(phi1,[3*n_nodos,1]);
dphi0 = reshape(dphi0,[3*n_nodos,1]);
dphi1 = reshape(dphi1,[3*n_nodos,1]);
%parte explicita
cte = (4/9)*(dt^2/(rho));
Fexp = Fg + K*(phi0 - 2*phi1) + D*(dphi0 - 2*dphi1);
phi = -(1/3)*phi0 + (4/3)*phi1 - (2/9)*dt*dphi0 + (8/9)*dt*dphi1 + cte*Minv*Fexp;
%restriccion
[Cphi,J] = fun_C(phi,C,Mcons,A_b,n_nodos,n_conds);
C0((end-(length(u)-1)):end) = u; den = C0;  den(abs(den) < 10^-6) = Inf;
error = Cphi-C0; error_max = 2; n_iter = 0;
while error_max > 0.5 || n_iter < 1
    MinvJt = cte*Minv*J';
    A = round(J*MinvJt,12);
    dlt_lambda = A\error;
    dlt_phi = -MinvJt*dlt_lambda;
    phi = phi + dlt_phi;
    [Cphi,J] = fun_C(phi,C,Mcons,A_b,n_nodos,n_conds);
    error = Cphi-C0;
    error_max = 100*max(abs(error./den));
    n_iter = n_iter+1;
end
%actualizamos la velocidad
dphi = (0.5/dt)*(3*phi - 4*phi1 + phi0);

