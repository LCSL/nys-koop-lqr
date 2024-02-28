% Code accompanying the IEEE TAC Paper "Optimal construction of Koopman eigenfunctions for prediction and control"
% by Milan Korda and Igor Mezic 

% Milan Korda, March 2020

clear all
close all

addpath('./Functions');
addpath('./Functions/qpOASES-3.1.0/interfaces/matlab')
try
    [~] = qpOASES(1,0,-1,1);
    QPsolver = 'qpoases';
catch        
    warning('qpOASES not installed properly, using quadprog as a backup')
    warning('To install it qpOASES, run ./Functions/qpOASES-3.1.0/interfaces/matlab/make.m')
    warning('If this does not work, see https://www.coin-or.org/qpOASES/doc/3.0/manual.pdf')
    QPsolver = 'quadprog';
end



%% Parameters
% Set to 0 if you want to reoptimize the eigenvalues
USE_PRECOMPUTED_LAMBDAS = 1; 

% Total number of eigenfunctions use
%Ns_efun = Ns_efun(1:6)
Ns_efun = int64(logspace(1, 2.3, 20)); 

all_rmses = [];
seed = 2141442;
rng(seed)
%% Damped Duffing oscillator dynamics
delta = 0.5;
beta = -1;
alpha = 1;
f = @(t,x)( 0.5 * [2*x(2,:) ; -delta*2*x(2,:) - 2*x(1,:)*(beta+alpha*(2*x(1,:)).^2)] ); % Uncontrolled
f_u = @(t,x,u)( 0.5 * [2*x(2,:) ; -delta*2*x(2,:) - 2*x(1,:).*(beta+alpha*(2*x(1,:)).^2) + u]  ); % Controlled

n = 2; % Number of states
m = 1; % Number of control inputs


%% Discretize

% Sampling period
deltaT = 0.01;

%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );


%% Collect trajectories
trajLen = 500; % Lenght of tracetories
Ntraj = 100; % Number of trajectories
Traj = cell(1,Ntraj); % Cell array of trajectories
for j = 1:Ntraj
    textwaitbar(j, Ntraj, "Collecting data without control")
    
    xx = randn(n,1);
    xx = xx / norm(xx); % Unitial conditions on unit circle
    for i = 1:trajLen-1
        xx = [xx f_ud(0,xx(:,end),0)];
    end
    Traj{j} = xx;
end

% Vectorize data
F_vec = cell(1,n);
for i = 1:n
    F_vec{i} = [];
    for j = 1:numel(Traj)
        F_vec{i} = [F_vec{i} ; Traj{j}(i,:).'];
    end
end

%% DMD
X_DMD = []; Y_DMD = [];
for i = 1:numel(Traj)
    X_DMD = [X_DMD, Traj{i}(:,1:end-1)];
    Y_DMD = [Y_DMD, Traj{i}(:,2:end)];
end
A_DMD = Y_DMD * pinv(X_DMD);
% writematrix(X_DMD, "duffing_x_unforced.csv")
% writematrix(Y_DMD, "duffing_y_unforced.csv")

%% Collect dta with control
disp('Collecting with control')
Nsim = 200;

% Random forcing
Ubig = 2*rand([Nsim Ntraj]) - 1;

% Generate trajectories with control
tic
% Initial conditions
Ntraj = numel(Traj);
Xcurrent = [];
for i = 1:numel(Traj)
    Xcurrent = [Xcurrent, Traj{i}(:,1) ];
end
X_Data = zeros(n,Nsim+1,Ntraj);
U_Data = zeros(m,Nsim,Ntraj);
X_Data(:,1,:) = Xcurrent;
X = []; Y = []; U = [];
for i = 1:Nsim
    textwaitbar(i,Nsim,"Collecting data with control")
    Xnext = f_ud(0,Xcurrent,Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    
    U_Data(:,i,:) = Ubig(i,:); % Record data trajectory - by trajectory
    X_Data(:,i+1,:) = Xnext;
    Xcurrent = Xnext;
end
% writematrix(X, "duffing_x_forced.csv")
% writematrix(U, "duffing_u_forced.csv")
% writematrix(Y, "duffing_y_forced.csv")

inds = false(1,Ntraj);
for i = 1:Ntraj
    if( ~any(any(abs(X_Data(:,:,i)) > 3)) && ~any(any(isnan(X_Data(:,:,i)))) )
        inds(i) = 1;
    end
end
% Convert to cell array
Xtraj = cell(1,Ntraj); Utraj = cell(1,Ntraj);
for i = 1:Ntraj
    Xtraj{i} = X_Data(:,:,i);
    Utraj{i} = U_Data(:,:,i);
end
Xtraj = Xtraj(inds); Utraj = Utraj(inds); Ntraj = numel(Xtraj);
validate_forecasts = false;
if ~validate_forecasts
    Ns_efun = 90;
end
for N_efun = Ns_efun                                                                    
    %% Eigenvalues

    % if(USE_PRECOMPUTED_LAMBDAS)
    %     load('lam_opt_duffing.mat','LAM_OPT')
    % else
    % DMD eigenvalues
    lam_dt = eig(A_DMD); 

    % Generate lattice of from DMD eigenvalues
    deg = 20;%4
    pows = monpowers(numel(lam_dt),deg);    
    lam_dt_lattice = [];
    for i = 1:size(pows,1)
        lam_dt_lattice = [ lam_dt_lattice ;  prod(lam_dt.'.^pows(i,:)) ];
    end   

    LAM_OPT{1} = lam_dt_lattice(1:N_efun/2);
    LAM_OPT{2} = lam_dt_lattice(1:N_efun/2);

    % end

    %% Optimization of Lambdas

    tic
    if(USE_PRECOMPUTED_LAMBDAS == 0)

        % Initial lambda
        lam0 = [real(LAM_OPT{1}); imag(LAM_OPT{1})];


        % Optiization function
        optimize_lam_fun_1 = @(x)(eigOptim_grad(x,F_vec{1},Ntraj,trajLen));
        optimize_lam_fun_2 = @(x)(eigOptim_grad(x,F_vec{2},Ntraj,trajLen));

        % Options
        options = optimoptions('fmincon','Display','iter','MaxFunEvals',1e6, 'SpecifyObjectiveGradient',true);

        % Optimize
        lam_opt_1 = fmincon(optimize_lam_fun_1,lam0,[],[],[],[],-1.05*ones(numel(lam0),1),1.05*ones(numel(lam0),1),[],options);
        lam_opt_2 = fmincon(optimize_lam_fun_2,lam0,[],[],[],[],-1.05*ones(numel(lam0),1),1.05*ones(numel(lam0),1),[],options);

        % Get back minimizers
        lam_opt_realImag{1} = lam_opt_1;
        lam_opt_realImag{2} = lam_opt_2;
        LAM_OPT{1} = lam_opt_1(1:N_efun/n) + 1i*lam_opt_1((N_efun/n)+1:end);
        LAM_OPT{2} = lam_opt_2(1:N_efun/n) + 1i*lam_opt_2((N_efun/n)+1:end);

        % Comapre residuals
        for i = 1:n
            residual_optim(i) = 100 * sqrt(eigOptim_grad(lam_opt_realImag{i},F_vec{i},Ntraj,trajLen)) / norm(F_vec{i}) ;
            fprintf('Regression residual (AFTER lambda optimization) for x_%d = %f %% \n',i,  residual_optim(i))
        end
        fprintf('Optimization time = %f \n',toc);
    end



    %% Regression for the boundary functions g

    lam_powers_traj = cell(n,numel(LAM_OPT{1}));
    for i = 1:n
        for j = 1:numel(LAM_OPT{i})
            lam_powers_traj{i,j} = bsxfun(@power,LAM_OPT{i}(j),0:trajLen-1);
        end
    end

    % Powers of eigenvalues along the trajectory
    lam_powers_traj = cell(n,numel(LAM_OPT{1}));
    for i = 1:n
        for j = 1:numel(LAM_OPT{i})
            lam_powers_traj{i,j} = bsxfun(@power,LAM_OPT{i}(j),0:trajLen-1);
        end
    end


    %Build linear operators mapping g0 to phi
    Lbig = cell(1,n);
    Lbig = [];
    L = cell(n,numel(LAM_OPT{i}));
    for i = 1:n
        Lbig{i} = [];
        for j = 1:numel(LAM_OPT{i})
            L{i,j} = bdiag(lam_powers_traj{i,j}.',Ntraj);
            Lbig{i} = [Lbig{i}, L{i,j}];
        end
    end


    for i = 1:n
        % g00{i} = Lbig{i} \ F_vec{i}; % Concatenated initial values
        g00{i} = (Lbig{i}'*Lbig{i}) \ (Lbig{i}'*F_vec{i}); % Seems to be numerically more stable
        for j = 1:numel(LAM_OPT{i})
            g0{i,j} = g00{i}((j-1)*Ntraj+1:j*Ntraj); % Initial values for each eigenfunction separately
        end
    end



    for i = 1:n
        residual_befre(i) = 100*norm(Lbig{i}*g00{i}-F_vec{i})/norm(F_vec{i});
        fprintf('Regression residual for x_%d = %f %% \n',i, residual_befre(i)  );
    end

    %% Convert trajectopries to a single vector
    X = [];
    for j = 1:numel(Traj)
        X = [X Traj{j}];
    end

    %% Values of eigenfunctions

    Val_phi = cell(n,numel(LAM_OPT{1}));
    for i = 1:n
        for j = 1:numel(LAM_OPT{1})
            Val_phi{i,j} = (L{i,j}*g0{i,j}).';
        end
    end


    %% Interpolate

    VAL_PHI_CONCAT = [];
    for i = 1:size(Val_phi,1)
         textwaitbar(i, size(Val_phi,1), "Interpolating")

        for k = 1:size(Val_phi,2)
            phi{i,k} = scatteredInterpolant(X',Val_phi{i,k}');
            phi{i,k} = @(x)(phi{i,k}(x')');
            VAL_PHI_CONCAT = [VAL_PHI_CONCAT ; Val_phi{i,k} ];
        end
    end

    % Vectorize
    ind = 1;
    phi_vec = cell(1,N_efun);
    for i = 1:n
        for j = 1:numel(LAM_OPT{i})
            phi_vec{ind} = phi{i,j};
            ind = ind + 1;
        end
    end


    %% Lifting function
    disp('Computing projection on the span of eigenfunctions')
    liftFun = @(x)( liftFun_function(x,phi_vec) );

    %% Matrices A and C
    C = bdiag(ones(1,N_efun/n),n);  % Could be obtained through regression as well
    A = diag(reshape(cell2mat(LAM_OPT),[],1)); % Put all eigenvalues on the diagonal


    %% Regression for B
    % Setup regression problem min ||Q * vec(Blift) - b||_2^2 corresponding to
    % min_{Blift} \sum_{j=1}^{Nsim} ||xpred - xtrue||^2 where
    % xpred = C*A^j*x0_lift + sum_{k=0}^{j-1} kron(u_k',C*A^{j-k-1})*vec(BLift)
    % (And this sum over all trajectories)
    Obs =  obsvk(sparse(A),C,Nsim+1); % Not an efficient implementation
    b = [];
    nc = size(C,1);
    Q = zeros(Ntraj*Nsim*nc,size(A,1)*m);
    for q = 1 : Ntraj
        textwaitbar(q,Ntraj,"Building regression matrces for B")

        x0 = Xtraj{q}(:,1);
        Obsx0 = Obs*liftFun(x0);
        for j = 1:Nsim
            b = [b ; Obsx0( j*nc + 1 : (j+1)*nc, : ) - Xtraj{q}(:,j+1)] ;
            tmp = 0;
            for k = 0 : j-1
                kprime = j - k -1;
                tmp = tmp + kron(Utraj{q}(:,k+1)',Obs(kprime*nc + 1 : (kprime+1)*nc,:));
            end
            Q((q-1)*Nsim*nc + (j-1)*nc + 1 : (q-1)*Nsim*nc + j*nc,:) = tmp;
        end
    end
    b = -b;

    B = pinv(Q)*b; % = vec(Blift)
    B = reshape(B,size(A,1),m); % Unvectorize
%     writematrix(full(A), strcat('A_', int2str(N_efun), '.csv'))
%     writematrix(full(B), strcat('B_', int2str(N_efun), '.csv'))
%     writematrix(full(C), strcat('C_', int2str(N_efun), '.csv'))


    %% Prediction comparison with control
    close all
    Tpred = 2;
    Npred = Tpred / deltaT;
    rmses = [];
    if validate_forecasts
        n_seeds = 200;
    else
        n_seeds = 1;
    end
    for seed = (1:n_seeds)
        rng(seed)
        fprintf("Curr seed and m: "+ seed + " " + N_efun + "\n")
        x0 = rand_sphere(n,1);
        % x0 = [-0.9685; -0.2689];
        Xlift = liftFun(x0);

        Xtrue = x0;
        u_dt = @(i)((-1).^(round(i/30)));

        for i = 1:Npred
            Xlift = [Xlift, A*Xlift(:,end) + B*u_dt(i)];
            Xtrue = [Xtrue, f_ud(0,Xtrue(:,end),u_dt(i))];
        end
        Xpred = real(C*Xlift);
        curr_rmse = sqrt(sum((Xpred - Xtrue).^2, "all") / sum(Xtrue.^2, "all")) * 100;
        rmses = [rmses;
                 curr_rmse];

    end
    all_rmses = [all_rmses, rmses];
end
if validate_forecasts
    writematrix(all_rmses, "all_rmses_eigfuns.csv");
end
R_lqr = eye(size(B, 2), size(B, 2));
Q_lqr = full(C' * C);

[K_lqr, ~, ~] = dlqr(A, B, Q_lqr, R_lqr);
x0 = [0.9; 0.8];
reference = [0; 0];
Xlift = liftFun(x0);
ref_lift = liftFun(reference);
Xtrue = x0;
u_dt = @(i)((-1).^(round(i/30)));
us = [];
for i = 1:5 * Npred
    u_dt = K_lqr * (ref_lift - Xlift);
    disp(u_dt)
    us = [us, u_dt];
    Xnew = f_ud(0,Xtrue(:,end),real(u_dt));
    Xtrue = [Xtrue, Xnew];
    Xlift = liftFun(Xnew);
end
scatter(0, 0, 'filled');
hold on
rec_r = real(C * ref_lift);
scatter(rec_r(1), rec_r(2), 'MarkerEdgeColor', [0 0 .5]);

for i = (1:size(Xtrue, 2))
    scatter(Xtrue(1, i), Xtrue(2, i), 'MarkerEdgeColor', [0 .5 .5]);
    hold on
    pause(0.00001)
end
hold off
% 
% 
% figure
% lw = 5;
% plot([0:Npred]*deltaT,Xtrue(1,:),'linewidth',lw); hold on
% plot([0:Npred]*deltaT,Xpred(1,:),'--r','linewidth',lw);
% LEG = legend('True','Predicted'); set(LEG,'fontsize',25,'interpreter','latex','location','southeast')
% ylim([min(min(Xtrue(1,:)), min(Xpred(1,:)))-0.1,max(max(Xtrue(1,:)), max(Xpred(1,:)))+0.1])
% set(gca,'FontSize',20);
% xlabel('Time [s]','interpreter','latex','fontsize',30);
% ylabel('$x_1$','interpreter','latex','fontsize',35);
% 
% 
% figure
% plot([0:Npred]*deltaT,Xtrue(2,:),'linewidth',lw); hold on
% plot([0:Npred]*deltaT,Xpred(2,:),'--r','linewidth',lw);
% LEG = legend('True','Predicted'); set(LEG,'fontsize',25,'interpreter','latex','location','southeast')
% set(gca,'FontSize',20);
% xlabel('Time [s]','interpreter','latex','fontsize',30);
% ylabel('$x_2$','interpreter','latex','fontsize',35);
% 
% disp('Press any key for feedback control')
% pause()
% 
% %% Feedback control using Koopman MPC 
% 
% Qy = diag([1;0.1]);
% R = 0.0001;
% Tpred = 1;
% Np = round(Tpred / deltaT);
% xlift_min = nan(N_efun,1);
% xlift_max = nan(N_efun,1);
% 
% % Get the Koopman MPC controller
% koopmanMPC  = getMPC(A,B,C,0,Qy,R,Qy,Np,-1, 1, xlift_min, xlift_max,QPsolver);
% 
% 
% Tsim_CL = 15;
% Nsim_CL = Tsim_CL / deltaT;
% 
% % Define the reference
% yrr = zeros(2,Nsim_CL);
% yrr(:,1:Nsim_CL/3) = repmat([-0.5 ; 0],1,Nsim_CL/3);
% yrr(:,Nsim_CL/3+1:2*Nsim_CL/3) = repmat([0 ; 0],1,Nsim_CL/3);
% yrr(:,2*Nsim_CL/3+1:end) = repmat([0.25 ; 0],1,Nsim_CL/3);
% 
% % Intial condition for the closed-loop simulation
% x0 = [0.5;0];
% 
% XX_koop = x0;
% x_koop = x0;
% 
% UU_koop = [];
% 
% % Closed-loop simultion start
% tic
% for i = 1:Nsim_CL
% 
%     textwaitbar(i,Nsim_CL,"Closed-loop simulation")
% 
%     yr = yrr(:,i);
%     % Lift
%     xlift = liftFun(x_koop);
% 
%     % Get control
%     U_koop = koopmanMPC(xlift,yr);
% 
%     % CL dynamics
%     x_koop = f_ud(0,x_koop,U_koop(1));
% 
%     % Store
%     XX_koop = [XX_koop x_koop];
%     UU_koop = [UU_koop U_koop(1)];
% end
% 
% % Plot
% close all
% lw = 5;
% figure
% plot([0:Nsim_CL-1]*deltaT,yrr(1,:),'--r','linewidth',lw); hold on
% plot([0:Nsim_CL]*deltaT,XX_koop(1,:),'linewidth',lw,'color',[0 0 0.8]);
% set(gca,'fontsize',20)
% xlabel('Time [s]','interpreter','latex','fontsize',30);
% ylabel('$x_1$','interpreter','latex','fontsize',30);
% LEG = legend('Reference','State $x_1$'); set(LEG,'interpreter','latex','fontsize',30,'location','southeast')
% 
% 
% figure
% plot([0:Nsim_CL]*deltaT,XX_koop(2,:),'linewidth',lw,'color',[0 0 0.8]); hold on
% set(gca,'fontsize',20)
% xlabel('Time [s]','interpreter','latex','fontsize',30);
% ylabel('$x_2$','interpreter','latex','fontsize',30);
% LEG = legend('State $x_2$'); set(LEG,'interpreter','latex','fontsize',30,'location','southeast')
% 
% 
% figure
% plot(XX_koop(1,1:Nsim_CL/3),XX_koop(2,1:Nsim_CL/3),'color',[0 0 0.8],'linewidth',lw); hold on
% plot(XX_koop(1,Nsim_CL/3+1:2*Nsim_CL/3),XX_koop(2,Nsim_CL/3+1:2*Nsim_CL/3),'color',[0.8 0 0],'linewidth',lw); hold on
% plot(XX_koop(1,2*Nsim_CL/3+1:end),XX_koop(2,2*Nsim_CL/3+1:end),'color',[0 0.7 0],'linewidth',lw); hold on
% scatter([-0.5 0 0.25],[0 0 0],20*[10 10 10],[0 0 0],'filled')
% scatter(0.5,0,20*[10],[0.8 0 0],'filled')
% set(gca,'fontsize',20)
% xlabel('$x_1$','interpreter','latex','fontsize',30);
% ylabel('$x_2$','interpreter','latex','fontsize',30);
% 
% 
% figure
% h1 = plot([0:Nsim_CL-1]*deltaT,-1*ones(1,Nsim_CL),'--k','linewidth',3); hold on
% h2 = plot([0:Nsim_CL-1]*deltaT,1*ones(1,Nsim_CL),'--k','linewidth',3);
% h3 = plot([0:Nsim_CL-1]*deltaT,UU_koop,'linewidth',lw,'color',[0 0 0.8]);
% set(gca,'fontsize',20)
% xlabel('Time [s]','interpreter','latex','fontsize',30);
% ylabel('$u$','interpreter','latex','fontsize',30);
% ylim([-1.1,1.1])
% LEG = legend([h3,h1],{'Control input','Constraints'}); set(LEG,'interpreter','latex','fontsize',30,'location','southeast')
