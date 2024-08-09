% Author: Edoardo Caldarelli, IRI, CSIC-UPC, Barcelona
% ecaldarelli@iri.upc.edu
% Based on code by Franco Coltraro, IRI, CSIC-UPC, Barcelona
% July 2024
close all
clear all

m = cast(100, 'int64');
%%%% Import Koopman python modules
pyrun('import numpy as np')
pyrun('import random')
pyrun('import pickle')
for seed = (0:49)
    rmses_nys = [];
    rmses_splines = [];
    errors_lqr_nys = [];
    errors_lqr_splines = [];

    for kapprox = ["nystrom", "splines"]
        disp(seed)
        pyrun("np.random.seed(" + seed + ")")
        pyrun("random.seed(" + seed + ")")
        
        pyrun("import importlib.util, sys")
        pyrun("spec = importlib.util.spec_from_file_location('regressors', '/home/ecaldarelli/PycharmProjects/koopman_mpc/regressors.py')")
        pyrun('module = importlib.util.module_from_spec(spec)')
        pyrun("sys.modules['regressors'] = module")
        pyrun("spec.loader.exec_module(module)")
        pyrun('from regressors import ThreeDimensionalKernel, KoopmanNystromRegressor, KoopmanRegressor')
        
        %pyrun("spec = importlib.util.spec_from_file_location('regressors_deprecated', '/Users/edoardocaldarelli/Documents/phd/koopman_mpc/regressors_deprecated.py')")
        pyrun('module = importlib.util.module_from_spec(spec)')
        pyrun("sys.modules['regressors_deprecated'] = module")
        pyrun("spec.loader.exec_module(module)")
        pyrun('from regressors_deprecated import KoopmanSplineRegressor')
        rng(0);
        
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
    
        delta_us = [];
        phiPositions = {};  phiPositions{1} = phi0;
        phiVelocities = {}; phiVelocities{1} = dphi0;
        %para iniciar el integrador de segundo orden
        delta_u = 0; %array(pyrun("res = arg1 @ (arg2 - arg3)", "res", arg1=K_lqr, arg2=lifted_ref, arg3=lifted_state));
        u(:,2) = u(:, 1) + delta_u;
        for i = (1:100)
            [phi1,dphi1] = explicit_grispun(dt,phi0,dphi0,...
                        Fg,rho,A_b,Mlum,Minv,Mcons,K,D,C,...
                        n_conds,n_nodos,Cphi0,u(:,2));
            phi0 = phi1;
            dphi0 = dphi1;
        end
        
        
                
        phiPositions{2} = reshape(phi1,[n_nodos,3]);
        phiVelocities{2} = reshape(dphi1,[n_nodos,3]);
        
        controlled_node1 = zeros(3,tf);
        tmp = reshape(phi0,[n_nodos,3]);
        controlled_node1(:,1) = tmp(controlados(1),:);
        
    
        ref = pyrun("res = np.loadtxt(f'/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/{arg1}/data/reference_lqr.csv').reshape([-1, 1])", "res", arg1=kapprox);
        tic
        % TODO careful with ordering of the nodes!!!!!!!!!!
        K_lqr = pyrun("res = np.loadtxt(f'/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/{arg1}/data/REBUTTAL_K_lqr_seed_" + seed + "_m_" + m + ".csv')", "res", arg1=kapprox);
        fid = py.open("/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/" + kapprox + "/data/REBUTTAL_regressor_seed_" + seed + "_m_" + m + ".npy", 'rb');
        regressor = py.pickle.load(fid);
        lifted_ref = pyrun("res = arg1.lift(arg2)", "res", arg1=regressor, arg2=ref);
        phi0_reshaped = reshape(reshape(phi0, n_nodos, 3)', n_nodos * 3, 1); 
        phi0_np = pyrun("res = np.array(arg1).reshape([-1, 1])", "res", arg1=phi0_reshaped);
        lifted_state = pyrun("res = arg1.lift(arg2)", "res", arg1=regressor, arg2=phi0_np);
        phiini = phi0;
        error_lqr = 0;
        for tt=1:(tf-1)  
            disp(tt)
            phi1_reshaped = reshape(reshape(phi1, n_nodos, 3)', n_nodos * 3, 1);
            phi_np = pyrun("res = np.array(arg1).reshape([-1, 1])", "res", arg1=phi1_reshaped);
            lifted_state = pyrun("res = arg1.lift(arg2)", "res", arg1=regressor, arg2=phi_np);
            curr_rmse = sqrt(sum(double(pyrun("res = (arg2 - arg3)", "res", arg2=ref, arg3=phi_np)).^2));
            if kapprox == "nystrom"
                rmses_nys = [rmses_nys, curr_rmse];
            else
                rmses_splines = [rmses_splines, curr_rmse];
            end
            delta_u = double(pyrun("res = arg1 @ (arg2 - arg3)", "res", arg1=K_lqr, arg2=lifted_ref, arg3=lifted_state));
            delta_us = [delta_us, delta_u];
            error_lqr = error_lqr + 0.0075 * (phi1 - double(ref))' * (phi1 - double(ref)) + delta_u' * delta_u;
           
            u(:, tt + 2) = u(:, tt + 1) + delta_u;
            %resolvemos las EDO
            [phi,dphi] = imex_franco(dt,phi0,dphi0,phi1,dphi1,...
                                   Fg,rho,A_b,Minv,Mcons,K,D,C,...
                                   n_conds,n_nodos,Cphi0,u(:, tt+2));
            tmp = reshape(phi,[n_nodos,3]);
            controlled_node1(:,tt+1) = tmp(controlados(1),:);
    
            phiPositions{tt+2} = reshape(phi,[n_nodos,3]);
            phiVelocities{tt+2} = reshape(dphi,[n_nodos,3]);   
            %actualizamos
            phi0 = phi1; dphi0 = dphi1;
            phi1 = phi;  dphi1 = dphi;
        end
        % plot(rmses_nys)
        % plot(rmses_splines)
        disp(error_lqr)
        if kapprox == "nystrom"
            errors_lqr_nys = [errors_lqr_nys, error_lqr];
        else
            errors_lqr_splines = [errors_lqr_splines, error_lqr];
        end
        toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%POSTPROCESO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %pelicula
        % peli = ''; vel_vid = 1; 
        % hacerPeli(phiPositions,T,nodos_borde,peli,vel_vid)
        % close all
        writematrix(delta_us, "/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/REBUTTAL_delta_us"+ kapprox + "_seed_" + seed + ".csv")

    end
    
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
    writematrix(errors_lqr_nys, "/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/REBUTTAL_lqrcost_nys_seed_" + seed + ".csv")
    writematrix(errors_lqr_splines, "/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/REBUTTAL_lqrcost_splines_seed_" + seed + ".csv")
    
    writematrix(errors_lqr_splines, "/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/REBUTTAL_lqrcost_splines_seed_" + seed + ".csv")

    writematrix(rmses_nys, "/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/REBUTTAL_rmse_nys_seed_" + seed + ".csv")
    writematrix(rmses_splines, "/home/ecaldarelli/PycharmProjects/koopman_mpc/8x8_cloth_swing_xyz/sim_results/REBUTTAL_rmse_splines_seed_" + seed + ".csv")
end