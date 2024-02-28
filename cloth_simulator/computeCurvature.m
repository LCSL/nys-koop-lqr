function out = computeCurvature(X,T,theReferenceElement,Minv)

%elemento
w = theReferenceElement.IPweights;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;
Nxieta=theReferenceElement.Nxieta;
nodos_elemento = size(theReferenceElement.nodesCoord,1);

%mallado
n_elementos = size(T,1);
n_nodos = size(X,1);
%creamos las matrices del sistema
Kgauss = zeros([n_nodos,1]);

for i=1:n_elementos
    Te=T(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element   
    Ke = zeros([nodos_elemento,1]); 
    for k=1:length(w)
        %calculo de las derivadas en el elemento
        N_k = N(k,:);
        Nxi_k = Nxi(k,:); Neta_k =  Neta(k,:); Nxieta_k =  Nxieta(k,:);
        phi_xi = (Nxi_k*Xe); phi_eta = (Neta_k*Xe);  phi_xieta = (Nxieta_k*Xe);
        %normal
        nu = cross(phi_xi,phi_eta); nu = nu/norm(nu);
        %forma fundamental
        E = (phi_xi)*(phi_xi)'; G = (phi_eta)*(phi_eta)'; F = (phi_xi)*(phi_eta)';
        %segunda forma fundamental
        f = (phi_xieta)*(nu');
        %elemento de area
        dS = sqrt(abs(E*G - F^2))*w(k);
        %curvatura
        Ke = Ke + ((-f^2)/(E*G - F^2))*dS*N_k';
    end    
    Kgauss(Te) = Ke + Kgauss(Te);
end


out = Minv(1:n_nodos,1:n_nodos)*Kgauss;