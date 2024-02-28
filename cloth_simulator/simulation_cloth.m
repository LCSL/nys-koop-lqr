% Simulador de telas inextensibles 
% Franco Coltraro, IRI, CSIC-UPC, Barcelona
% 2019/10/8
close all
% clear

%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETROS SIMULACION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tiempo
H = 5;
dt = 0.05;
tf = H/dt+2;
times = dt*(0:tf);
%parametros de la superficie
rho = 0.3; %densidad de masa inercial
delta = 0.3; %masa gravitatoria (aerodinamica)
theta = 0.0005; %resistencia  doblarse: rigidez
alfa = 0.6;  %amortiguamiento oscilaciones lentas
beta = 0.01*theta; %amortiguamiento oscilaciones rapidas
c = 3000; %cizallamiento (c=0 es inextensibilidad)
%parametros del objeto del choque
eps0 = 0.002; %distancia al objeto (i.e. grosor tela)
mu = [0.5,0,0.75]; %friccion suelo, objeto extra y la tela consigo misma
%parametros no fisicos
% rho = 0.2; %densidad de masa
% theta = 0; %resistencia  doblarse: rigidez;
% alfa = 1.5;  %amortiguamiento oscilaciones lentas (~rozamiento con el aire)
% beta = 0; %amortiguamiento oscilaciones rapidas
% % rho = 0.1; %densidad de masa
% % theta = 10; %resistencia  doblarse: rigidez;
% % alfa = 0.15;  %amortiguamiento oscilaciones lentas (~rozamiento con el aire)
% % beta = 4; %amortiguamiento oscilaciones rapidas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MALLADO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%se crea a partir de una parametrizacion de la forma (f1(x,y),f2(x,y),f3(x,y))
%donde (x,y) estan en el rectangulo [ax,bx]x[ay,by]
ax = -0.5; bx = 0.5; ay = -0.5; by = 0.5; 
f1 = @(x,y) x; f2 = @(x,y) 0.75*sqrt(1 - x.^2) +0.1; f3 = @(x,y) y;
npx = 8; npy = 8; n_nodos = npx*npy; %numero de puntos
[X,T] = CreateMesh([ax,bx,ay,by],npx,npy,f1,f2,f3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOUNDARY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%deteccion de los puntos del borde para poder imponer isometria tambien
%alli
[Xb,Tb,nodos_borde] = GetBoundary(X,[ax,bx,ay,by],npx,npy,f1,f2,f3); 
n_nodos_bnd = size(Xb.Xb,1);
%interior nodes
nodes_int = setdiff(1:n_nodos,nodos_borde.nodes_bnd); 
%coordenadas de las esquinas
esquinas = [1, npx, npx*(npy-1)+1, npx*npy]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%VARIABLES DE ESTADO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Superficie inicial: ver CreateMesh
phi0 = X; 
%Velocidad inicial
dphi0 = [zeros([n_nodos,1]),zeros([n_nodos,1]),zeros([n_nodos,1])]; 
%%%%%%%%%%%%%%%%%%%%%%%%%VARIABLES DE ESTADO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%NODOS A CONTROLAR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%damos las coordenadas de los nodos de los que vamos a fijar su trayectoria
controlados = [npx*(npy-1)+1, npx*npy]; %dos esquinas
coord_controlados = [controlados, controlados+n_nodos, controlados+(2*n_nodos)];
%matriz para imponer las condiciones en el contorno
A_b = spalloc(length(coord_controlados),3*n_nodos,length(coord_controlados));
A_b(:,coord_controlados) = eye(length(coord_controlados));

%posicion de los nodos para el control
%rng('default')
u = zeros(length(coord_controlados),tf+1);
u(:,1) = reshape(X(controlados,:),[length(coord_controlados),1]);
delta_u = 0;

% select random frequencies
min_freq = 0.3;
max_freq = 0.6;
freqX = min_freq + (max_freq-min_freq)*rand; % [Hz]
freqY = min_freq + (max_freq-min_freq)*rand; % [Hz]
freqZ = min_freq + (max_freq-min_freq)*rand; % [Hz]

% select random phases
min_phase = -0.1;
max_phase = 0.1;
phaseX = min_phase + (max_phase-min_phase)*rand;
phaseY = min_phase + (max_phase-min_phase)*rand;
phaseZ = min_phase + (max_phase-min_phase)*rand;

% select random direction in Y-Z
min_angle_u = -60; % [deg]
max_angle_u = 60; % [deg]
angle_u = min_angle_u + (max_angle_u-min_angle_u)*rand; % [deg]

% select random amplitudes
min_amplX = 0.0;%0.04;
max_amplX = 0.0;%0.08;
amplX = min_amplX + (max_amplX-min_amplX)*rand;
amplY = -0.1 * cos(angle_u*pi/180);
amplZ = 0.1 * sin(angle_u*pi/180);
for tt=2:(tf+1)
        v = [0;0;-1;-1;0;0]; %Y
        delta_u = amplY*cos(2*pi*freqY*times(tt)+phaseY)*v + 0.0*randn(1)*v;
        v = [-1;-1;0;0;0;0]; % X
        delta_u = delta_u  + amplX*cos(2*pi*freqX*times(tt)+phaseX)*v + 0.00*randn(1)*v;
        v = [0;0;0;0;-1;-1]; % Z
        delta_u = delta_u + amplZ*cos(2*pi*freqZ*times(tt)+phaseZ)*v + 0.0*randn(1)*v; 
        u(:,tt) = u(:,tt-1) + delta_u;
end

% nodos1 = T(1,:);
% nodos2 = T(T(:,2)==npx,:);
% %matriz para imponer las condiciones en el contorno
% A_b = spalloc(6,3*n_nodos,24);
% unos = 0.25*ones([1 length(nodos1)]);
% A_b(1:3,[nodos1,nodos1+n_nodos,nodos1+(2*n_nodos)]) = blkdiag(unos,unos,unos);
% A_b(4:6,[nodos2,nodos2+n_nodos,nodos2+(2*n_nodos)]) = blkdiag(unos,unos,unos);
% posicion de los nodos para el control
% u = zeros(6,tf+1);
% u(:,1) = A_b*X(:);
% v = [0;-1;0;0;-1;0]; freq = 1;
% for tt=2:(tf+1)
%     if times(tt) < 4
%         u(:,tt) = u(:,tt-1);
%     elseif times(tt) >= 4 && times(tt) <= 6
%         u(:,tt) = u(:,tt-1) + 0.015*cos(2*pi*freq*times(tt)+0.1)*v;% + 0.008*randn([6,1]);
%     else
%         u(:,tt) = u(:,tt-1);
%     end 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%ELEMENTOS DE REFERENCIA%%%%%%%%%%%%%%%%%%%%%%%%%
%cuadrilateros bilineales
theReferenceElement = createReferenceElement();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%SYSTEM OF EQUATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculo de las matrices del sistema: masa Mlum Minv, rigidez K, 
%amortigualmiento D e inextensibilidad C (3-tensor)
[C,n_conds,Mlum,Minv,Mcons,D,K] = computeMatrices(X,T,...
                               nodes_int,nodos_borde,esquinas,...
                               [theta,alfa,beta],...
                               theReferenceElement);
[Cphi0,~] = fun_C(X,C,Mcons,A_b,n_nodos,n_conds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
%%%%%%%%%%%%%%%%%%%%GRAVEDAD Y FUERZAS EXTERNAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fg = Mlum*reshape([zeros([n_nodos,1]),...
                   zeros([n_nodos,1]),...
          -9.8*rho*ones([n_nodos,1])],[3*n_nodos,1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%SUPER BUCLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%almacenamiento de las trajectorias
phiPositions = {};  phiPositions{1} = phi0;
phiVelocities = {}; phiVelocities{1} = dphi0;

%para iniciar el integrador de segundo orden
[phi1,dphi1] = explicit_grispun(dt,phi0,dphi0,...
                Fg,rho,A_b,Mlum,Minv,Mcons,K,D,C,...
                n_conds,n_nodos,Cphi0,u(:,2));
            
phiPositions{2} = reshape(phi1,[n_nodos,3]);
phiVelocities{2} = reshape(dphi1,[n_nodos,3]);

tic
controlled_node1 = zeros(3,tf);
tmp = reshape(phi0,[n_nodos,3]);
controlled_node1(:,1) = tmp(controlados(1),:);
for tt=1:(tf-1)  
    disp(tt)
    %resolvemos las EDO
    [phi,dphi] = imex_franco(dt,phi0,dphi0,phi1,dphi1,...
                           Fg,rho,A_b,Minv,Mcons,K,D,C,...
                           n_conds,n_nodos,Cphi0,u(:,tt+2));
                       
    tmp = reshape(phi,[n_nodos,3]);
    controlled_node1(:,tt+1) = tmp(controlados(1),:);
             
    phiPositions{tt+2} = reshape(phi,[n_nodos,3]);
    phiVelocities{tt+2} = reshape(dphi,[n_nodos,3]);   
    %actualizamos
    phi0 = phi1; dphi0 = dphi1;
    phi1 = phi;  dphi1 = dphi;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%POSTPROCESO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pelicula
peli = ''; vel_vid = 1; 
%hacerPeli(phiPositions,T,nodos_borde,peli,vel_vid)

% % %Area
% areaError = computeAreaError(phiPositions,T,theReferenceElement);
% figure(2)
% plot(times,areaError)
% hold on
% mu = mean(areaError);
% hline = refline([0 mu]);
% hline.Color = 'r';
% hold off
% xlabel('Time')
% ylabel('Error')
% title('Relative total area error')
% disp('Error area medio:')
% disp(mu)
% ylim([0 10])

% %aristas
% aristasError = computeEdgesError(phiPositions,Tb.Tb);
% figure(3)
% plot(times,aristasError)
% hold on
% mu = mean(aristasError);
% hline = refline([0 mu]);
% hline.Color = 'r';
% hold off
% xlabel('Time')
% ylabel('Error')
% title('Relative boundary length error')
% ylim([0 10])

% % %energia
% energias = computeEnergies(phiPositions,phiVelocities,rho,Mlum,Fg,K);
% figure(4)
% plot(times,energias)
% hold on
% mu = mean(energias);
% disp('Energia media:')
% disp(mu)
% hline = refline([0 mu]);
% hline.Color = 'r';
% xlabel('Time')
% ylabel('Energy')
% ylim([min(energias)-1 max(energias)+1])

% %guardamos info
% matrices = struct('Mlum',Mlum,'Minv',Minv,'K',(1/theta)*K,'C',C);
% save('manta.mat','phi','T','nodos_borde','matrices','theReferenceElement')

save('/Users/edoardocaldarelli/Documents/phd/MPC/SIMULADORforCGPDM/LOG/log_cloth.mat','phiPositions','phiVelocities','u','controlados','esquinas', 'angle_u', 'freqY', 'freqZ')