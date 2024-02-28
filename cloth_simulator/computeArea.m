function A = computeArea(X,T,theReferenceElement)

%elemento
w = theReferenceElement.IPweights;
Nxi=theReferenceElement.Nxi;
Neta=theReferenceElement.Neta;

%mallado
n_elementos = size(T,1);

%creamos las matrices del sistema
A = zeros([n_elementos,1]);

for i=1:n_elementos
    Te=T(i,:); %index of the nodes in the element
    Xe=X(Te,:); %coordinates of the nodes in the element   
    Ae = 0; 
    for k=1:length(w)
        %calculo de las derivadas en el elemento
        Nxi_k = Nxi(k,:); Neta_k =  Neta(k,:);
        %elemento de area
        E = (Nxi_k*Xe)*(Nxi_k*Xe)';
        G = (Neta_k*Xe)*(Neta_k*Xe)';
        F = (Nxi_k*Xe)*(Neta_k*Xe)';
        dS = sqrt(abs(E*G - F^2))*w(k);
        %lo juntamos todo
        Ae = Ae + dS; %2-tensor   
    end    
    A(i) = Ae;
end