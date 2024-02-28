function theReferenceElement = createReferenceElement()
%Creates a struct with the information of the reference element, for
%quad discretizations. Only supports bilinear quads.
typeElement='QUA';
nodesCoord = [-1 -1; 1 -1; 1 1; -1 1];
z = (1/sqrt(3))*[-1 -1; 1 -1; 1 1; -1 1];
w = [1 1 1 1];
xi = z(:,1); eta=z(:,2);
%funciones de forma y sus derivadas
N = (1/4)*[(1-xi).*(1-eta),(1+xi).*(1-eta),(1+xi).*(1+eta),(1-xi).*(1+eta)];
Nxi  = (1/4)*[ eta - 1 , 1 - eta, 1 + eta, -1 - eta ];
Neta = (1/4)*[xi - 1, -1 - xi, 1 + xi, 1 - xi];
Nxixi = 0*Nxi;
Netaeta = 0*Neta;
Nxieta = (1/4)*[ 0*eta + 1 , -1 - 0*eta, 1 + 0*eta, -1 - 0*eta ];
%struct
theReferenceElement = struct('IPweights',w,'IPcoord',z,'N',N,'Nxi',Nxi,'Neta',Neta,...
                           'Nxixi',Nxixi,'Netaeta',Netaeta,'Nxieta',Nxieta,...
                           'type',typeElement,'nodesCoord',nodesCoord);
end