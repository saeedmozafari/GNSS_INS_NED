
% Constants
deg_to_rad = 0.01745329252;
% rad_to_deg = 1/deg_to_rad;
% micro_g_to_meters_per_second_squared = 9.80665E-6;

% Position initialization error (m; N,E,D)
initialization_errors.delta_r_eb_n = [4;2;3];
% Velocity initialization error (m/s; N,E,D)
initialization_errors.delta_v_eb_n = [0.05;-0.05;0.1];
% Attitude initialization error (deg, converted to rad; @N,E,D)
initialization_errors.delta_eul_nb_n = [-0.05;0.04;1]*deg_to_rad; % rad