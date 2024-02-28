function fn_export_state_control(name,num)
    % Load LOG
    load('LOG/log_cloth.mat')
    [num_nodes, ~] = size(phiPositions{1});
    num_steps = length(phiPositions);

    input_dim = length(controlados);

    free_nodes = linspace(1,num_nodes,num_nodes);
    free_nodes(controlados) = [];


    % Build State and Input samples
    state = zeros(num_nodes*3, num_steps);
    input = zeros(input_dim*3, num_steps);
    for t=1:num_steps
        state(:,t) = reshape(phiPositions{t}', [num_nodes*3,1]);
        input(:,t) = reshape(phiPositions{t}(controlados,:)', [input_dim*3,1]);
    end
    state = state(:,1:end-1);
    input = input(:,2:end)-input(:,1:end-1);


    % Export DATA
    csvwrite(strcat('/Users/edoardocaldarelli/Documents/phd/MPC/SIMULADORforCGPDM/LOG/CONTROL_DATA/state_samples_',name,'_',num2str(num),'.csv'), state')
    csvwrite(strcat('/Users/edoardocaldarelli/Documents/phd/MPC/SIMULADORforCGPDM/LOG/CONTROL_DATA/input_samples_',name,'_',num2str(num),'.csv'), input')
    csvwrite(strcat('/Users/edoardocaldarelli/Documents/phd/MPC/SIMULADORforCGPDM/LOG/CONTROL_DATA/angle_',name,'_',num2str(num),'.csv'), angle_u)
    csvwrite(strcat('/Users/edoardocaldarelli/Documents/phd/MPC/SIMULADORforCGPDM/LOG/CONTROL_DATA/frequencies_',name,'_',num2str(num),'.csv'), [freqY; freqZ])
end
