close all; clc;
%% ========================================================================
Set_Initialization_Error;
%% ========================================================================
load('IMU_meas_60sec_car.mat');
%% ========================================================================
% Initialize true navigation solution
old_time = in_profile(1,1);
old_true_L_b = in_profile(1,2);
old_true_lambda_b = in_profile(1,3);
old_true_h_b = in_profile(1,4);
old_true_v_eb_n = in_profile(1,5:7)';
old_true_eul_nb = in_profile(1,8:10)';
old_true_C_b_n = Euler_to_CTM(old_true_eul_nb)';

% Initialize estimated navigation solution
[old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n] =...
    Initialize_NED(old_true_L_b,old_true_lambda_b,old_true_h_b,...
    old_true_v_eb_n,old_true_C_b_n,initialization_errors);
%% ========================================================================
out_errors(1,1) = old_time;
out_errors(1,2:4) = initialization_errors.delta_r_eb_n';
out_errors(1,5:7) = initialization_errors.delta_v_eb_n';
out_errors(1,8:10) = initialization_errors.delta_eul_nb_n';

% Progress bar
dots = '....................';
bars = '||||||||||||||||||||';
rewind = '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b';
fprintf(strcat('Processing: ',dots));
progress_mark = 0;
progress_epoch = 0;
%% ========================================================================
% Main loop
no_epochs = 6000;
for epoch = 2:no_epochs
    
    % Update progress bar
    if (epoch - progress_epoch) > (no_epochs/20)
        progress_mark = progress_mark + 1;
        progress_epoch = epoch;
        fprintf(strcat(rewind,bars(1:progress_mark),...
            dots(1:(20 - progress_mark))));
    end % if epoch
    
    % Input data from motion profile
    time = IMU_meas(epoch-1,1);
    
    % Time interval
    tor_i = time - old_time;
    
    meas_f_ib_b = IMU_meas(epoch-1,2:4)';
    meas_omega_ib_b = IMU_meas(epoch-1,5:7)';
    
    % Update estimated navigation solution
    [est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n] = ...
        Nav_equations_NED(tor_i,old_est_L_b,old_est_lambda_b,old_est_h_b,...
        old_est_v_eb_n,old_est_C_b_n,meas_f_ib_b,meas_omega_ib_b);
    
    % Generate output profile record
    out_profile(epoch,1) = time;
    out_profile(epoch,2) = est_L_b;
    out_profile(epoch,3) = est_lambda_b;
    out_profile(epoch,4) = est_h_b;
    out_profile(epoch,5:7) = est_v_eb_n';
    out_profile(epoch,8:10) = CTM_to_Euler(est_C_b_n')';
    
    % Reset old values
    old_est_L_b = est_L_b;
    old_est_lambda_b = est_lambda_b;
    old_est_h_b = est_h_b;
    old_est_v_eb_n = est_v_eb_n;
    old_est_C_b_n = est_C_b_n;
    
    %% ========================================================================
    % Determine errors and generate output record
    true_L_b = in_profile(epoch,2);
    true_lambda_b = in_profile(epoch,3);
    true_h_b = in_profile(epoch,4);
    true_v_eb_n = in_profile(epoch,5:7)';
    true_eul_nb = in_profile(epoch,8:10)';
    true_C_b_n = Euler_to_CTM(true_eul_nb)';
    
    [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
        est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n,true_L_b,...
        true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
    
    old_time = time;
    old_true_L_b = true_L_b;
    old_true_lambda_b = true_lambda_b;
    old_true_h_b = true_h_b;
    old_true_v_eb_n = true_v_eb_n;
    old_true_C_b_n = true_C_b_n;
    
    out_errors(epoch,1) = time;
    out_errors(epoch,2:4) = delta_r_eb_n';
    out_errors(epoch,5:7) = delta_v_eb_n';
    out_errors(epoch,8:10) = delta_eul_nb_n';
end %epoch
% Complete progress bar
fprintf(strcat(rewind,bars,'\n'));

% Plot the input motion profile and the errors (may not work in Octave).
close all;
Plot_profile(in_profile);
Plot_errors(out_errors);
