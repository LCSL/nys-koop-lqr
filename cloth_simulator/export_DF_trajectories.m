clear all
close all
clc
% let's indicate with P1 and P2 the two controlled corners
% and with F1 and F2 the two free corners.

% load('LOG/log_expl.mat')
load('LOG/log_ctrl.mat')
% u: [P1x, P2x, P1y, P2y, P1z, P2z]'

state_samples = zeros(30,length(phiPositions));

for t=1:length(state_samples)
    % --- t free corners' positions -----------------------------------
    state_samples(1,t) = phiPositions{t}(esquinas(1),1); % F1x(t)
    state_samples(2,t) = phiPositions{t}(esquinas(1),2); % F1y(t)
    state_samples(3,t) = phiPositions{t}(esquinas(1),3); % F1z(t)
    state_samples(4,t) = phiPositions{t}(esquinas(2),1); % F2x(t)
    state_samples(5,t) = phiPositions{t}(esquinas(2),2); % F2y(t)
    state_samples(6,t) = phiPositions{t}(esquinas(2),3); % F2z(t)
    if t>1 % --- t-1 free corners' positions --------------------------
      state_samples(7,t) = phiPositions{t-1}(esquinas(1),1); % F1x(t-1)
      state_samples(8,t) = phiPositions{t-1}(esquinas(1),2); % F1y(t-1)
      state_samples(9,t) = phiPositions{t-1}(esquinas(1),3); % F1z(t-1)
      state_samples(10,t) = phiPositions{t-1}(esquinas(2),1); % F2x(t-1)
      state_samples(11,t) = phiPositions{t-1}(esquinas(2),2); % F2y(t-1)
      state_samples(12,t) = phiPositions{t-1}(esquinas(2),3); % F2z(t-1)
    end 
    if t>2 % --- t-1 free corners' positions --------------------------
      state_samples(13,t) = phiPositions{t-2}(esquinas(1),1); % F1x(t-2)
      state_samples(14,t) = phiPositions{t-2}(esquinas(1),2); % F1y(t-2)
      state_samples(15,t) = phiPositions{t-2}(esquinas(1),3); % F1z(t-2)
      state_samples(16,t) = phiPositions{t-2}(esquinas(2),1); % F2x(t-2)
      state_samples(17,t) = phiPositions{t-2}(esquinas(2),2); % F2y(t-2)
      state_samples(18,t) = phiPositions{t-2}(esquinas(2),3); % F2z(t-2)
    end 
    % --- t pinched corners' positions -----------------------------------
    state_samples(19,t) = phiPositions{t}(controlados(1),1); % P1x(t)
    state_samples(20,t) = phiPositions{t}(controlados(1),2); % P1y(t)
    state_samples(21,t) = phiPositions{t}(controlados(1),3); % P1z(t)
    state_samples(22,t) = phiPositions{t}(controlados(2),1); % P2x(t)
    state_samples(23,t) = phiPositions{t}(controlados(2),2); % P2y(t)
    state_samples(24,t) = phiPositions{t}(controlados(2),3); % P2z(t)
    if t>1 % --- t-1 pinched corners' positions --------------------------
      state_samples(25,t) = phiPositions{t-1}(controlados(1),1); % P1x(t-1)
      state_samples(26,t) = phiPositions{t-1}(controlados(1),2); % P1y(t-1)
      state_samples(27,t) = phiPositions{t-1}(controlados(1),3); % P1z(t-1)
      state_samples(28,t) = phiPositions{t-1}(controlados(2),1); % P2x(t-1)
      state_samples(29,t) = phiPositions{t-1}(controlados(2),2); % P2y(t-1)
      state_samples(30,t) = phiPositions{t-1}(controlados(2),3); % P2z(t-1)
    end
    if t>2 % --- t-2 pinched corners' positions --------------------------
      state_samples(31,t) = phiPositions{t-2}(controlados(1),1); % P1x(t-2)
      state_samples(32,t) = phiPositions{t-2}(controlados(1),2); % P1y(t-2)
      state_samples(33,t) = phiPositions{t-2}(controlados(1),3); % P1z(t-2)
      state_samples(34,t) = phiPositions{t-2}(controlados(2),1); % P2x(t-2)
      state_samples(35,t) = phiPositions{t-2}(controlados(2),2); % P2y(t-2)
      state_samples(36,t) = phiPositions{t-2}(controlados(2),3); % P2z(t-2)
    end 
end
pinched_indeces = 19:24;
input_samples = state_samples(pinched_indeces,2:end)-state_samples(pinched_indeces,1:end-1);
state_samples = state_samples(:,1:end-1);

% Select samples
start_i = 3;
end_i = 92;
state_samples = state_samples(:,start_i:end_i);
input_samples = input_samples(:,start_i:end_i);

figure()
% plot free corners' trajectories
for m=0:2
    subplot(3,2,1)
    plot(state_samples(1+m*6,:))
    hold on
    title('$F^1$','interpreter','latex')
    ylabel('$x$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,3)
    plot(state_samples(2+m*6,:))
    hold on
    ylabel('$y$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,5)
    plot(state_samples(3+m*6,:))
    hold on
    ylabel('$z$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,2)
    plot(state_samples(4+m*6,:))
    hold on
    title('$F^2$','interpreter','latex')
    ylabel('$x$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,4)
    plot(state_samples(5+m*6,:))
    hold on
    ylabel('$y$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,6)
    plot(state_samples(6+m*6,:))
    hold on
    ylabel('$z$','interpreter','latex')
    grid on
    axis tight
end
figure()
% plot controlled corners' trajectories
for m=0:2
    subplot(3,2,1)
    plot(state_samples(19+m*6,:))
    hold on
    title('$P^1$','interpreter','latex')
    ylabel('$x$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,3)
    plot(state_samples(20+m*6,:))
    hold on
    ylabel('$y$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,5)
    plot(state_samples(21+m*6,:))
    hold on
    ylabel('$z$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,2)
    plot(state_samples(22+m*6,:))
    hold on
    title('$P^2$','interpreter','latex')
    ylabel('$x$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,4)
    plot(state_samples(23+m*6,:))
    hold on
    ylabel('$y$','interpreter','latex')
    grid on
    axis tight
    subplot(3,2,6)
    plot(state_samples(24+m*6,:))
    hold on
    ylabel('$z$','interpreter','latex')
    grid on
    axis tight
end

csvwrite('LOG/state_samples.csv',state_samples')
csvwrite('LOG/input_samples.csv',input_samples')