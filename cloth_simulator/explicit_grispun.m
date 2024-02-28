function [phi1,dphi1] = explicit_grispun(dt,phi,dphi,Fg,rho,A_b,Mlum,Minv,Mcons,K,D,C,n_conds,...
                    n_nodos,C0,u)
%hacemos un reshape
phi = reshape(phi,[3*n_nodos,1]);
dphi = reshape(dphi,[3*n_nodos,1]);
%paso sin restricciones

w0 = (Fg + ((rho/dt^2)*Mlum - K)*phi + ((rho/dt)*Mlum - D)*dphi);
phi0 = (dt^2/rho)*Minv*w0;
%restriccion
[Cphi0,J] = fun_C(phi0,C,Mcons,A_b,n_nodos,n_conds);
C0((end-(length(u)-1)):end) = u;
error = Cphi0-C0; n_iter = 0;
while max(abs(error)) > 10^-6 || n_iter < 1
    MinvJt = (1/rho)*Minv*J';
    dlt_lambda = ((dt^2)*J*MinvJt)\error;
    dlt_phi = -(dt^2)*MinvJt*dlt_lambda;
    phi1 = phi0 + dlt_phi;
    [Cphi1,J] = fun_C(phi1,C,Mcons,A_b,n_nodos,n_conds);
    error = Cphi1-C0;
    phi0 = phi1;
    n_iter = n_iter+1;
end
%actualizamos la velocidad
dphi1 = (phi1-phi)/dt;

