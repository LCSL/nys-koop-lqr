function [C,n_conds,Mlum,Minv,Mcons,D,K] = computeMatrices(X,T,...
                                         nodes_int,nodos_borde,esquinas,...
                                         params,theReferenceElement)
%Matrices del sistema: masa Mlum Minv, rigidez K, amortigualmiento D
%e inextensibilidad C (3-tensor aplanado)
theta = params(1);
alfa = params(2);
beta = params(3);
%elemento
w = theReferenceElement.IPweights;
N=theReferenceElement.N;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;
nodos_elemento = size(theReferenceElement.nodesCoord,1);

%mallado
n_nodos = size(X,1);
n_elementos = size(T,1);

%creamos las matrices del sistema
M = spalloc(n_nodos,n_nodos,9*n_nodos);
L = spalloc(n_nodos,n_nodos,9*n_nodos);
Cu = spalloc(n_nodos^2,n_nodos,27*n_nodos);
Cv = spalloc(n_nodos^2,n_nodos,27*n_nodos);
Cuv = spalloc(n_nodos^2,n_nodos,27*n_nodos);

%Loop in elements
for i=1:n_elementos
    Te=T(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element  
    %elemental matrices
    Me = zeros(nodos_elemento,nodos_elemento); 
    Le = zeros(nodos_elemento,nodos_elemento); 
    Cu_e = zeros(nodos_elemento^2,nodos_elemento);
    Cv_e = zeros(nodos_elemento^2,nodos_elemento);
    Cuv_e = zeros(nodos_elemento^2,nodos_elemento);
    %Bucle en puntos de integracion
    for k=1:length(w)
        %calculo de las derivadas en el elemento
        N_k = N(k,:); Nxi_k = Nxi(k,:); Neta_k =  Neta(k,:);
        %elemento de area
        E = (Nxi_k*Xe)*(Nxi_k*Xe)';
        G = (Neta_k*Xe)*(Neta_k*Xe)';
        F = (Nxi_k*Xe)*(Neta_k*Xe)';
        dS = sqrt(abs(E*G - F^2))*w(k);
        %lo juntamos todo
        Me = Me + kron(N_k',N_k)*dS; %2-tensor masa
        Le = Le + (Nxi_k'*Nxi_k + Neta_k'*Neta_k)*dS; %2-tensor laplaciano 
        %las tres condiciones para preservar la metrica
        Cu_e = Cu_e + kron(Nxi_k',kron(Nxi_k',N_k))*dS; %3-tensor 
        Cv_e = Cv_e + kron(Neta_k',kron(Neta_k',N_k))*dS; %3-tensor 
        Cuv_e = Cuv_e + 0.5*(kron(Nxi_k',kron(Neta_k',N_k))+kron(Neta_k',kron(Nxi_k',N_k)))*dS; %3-tensor 
    end    
    %assembly of elemental 2-tensors
    M(Te,Te) = M(Te,Te) + Me;
    L(Te,Te) = L(Te,Te) + Le;
    %3 tensor
    indices = combvectores(Te,Te); 
    ii = indices(1,:); jj = indices(2,:);
    TeTe = (jj - 1)*n_nodos + ii;
    Cu(TeTe,Te) = Cu(TeTe,Te) + Cu_e; 
    Cv(TeTe,Te) = Cv(TeTe,Te) + Cv_e; 
    Cuv(TeTe,Te) = Cuv(TeTe,Te) + Cuv_e; 
end

%solo nos quedamos con las condiciones que tengan sentido (en los bordes)
ind_u = [nodes_int,nodos_borde.nodes_bnd1(2:end)',nodos_borde.nodes_bnd3(1:(end-1))'];
ind_v = [nodes_int,nodos_borde.nodes_bnd2(2:end)',nodos_borde.nodes_bnd4(1:(end-1))'];
ind_uv = [nodes_int,esquinas];
C = [Cu(:,ind_u),Cv(:,ind_v),Cuv(:,ind_uv)]; 
n_conds = size(C,2);
%matrices
M_lumped = sum(M,1); 
Mlum = sparse(blkdiag(diag(M_lumped),diag(M_lumped),diag(M_lumped)));
Minv = diag(1./M_lumped);
Mcons = sparse(blkdiag(Minv(ind_u,ind_u),Minv(ind_v,ind_v),Minv(ind_uv,ind_uv)));
C = C*Mcons;
Kben = sparse(blkdiag(L'*Minv*L,L'*Minv*L,L'*Minv*L));
Minv = sparse(blkdiag(Minv,Minv,Minv));
D = alfa*Mlum + beta*Kben;
K = theta*Kben;



