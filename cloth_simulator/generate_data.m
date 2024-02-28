total_num = 50;
for num=0:total_num-1
    % run cloth simulation
    simulation_cloth
    
    % export data
    name = 'cloth_swing'; % choose name
    fn_export_state_control(name,num) % export csv files
end