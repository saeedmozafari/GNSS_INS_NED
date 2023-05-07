function [L_est,lambda_est,h_est,est_v_eb_n_new,est_C_b_n_new,est_IMU_bias_new,P_matrix_new] =...
    LC_KF_NED(GNSS_r_eb_n,GNSS_v_eb_n,tor_s,L,lambda,h,est_v_eb_n,est_C_b_n_old,est_IMU_bias_old,...
    P_matrix_old,meas_fBIB,meas_wBIB,LC_KF_config,lGBB)

D2R = 0.01745329252;%         convert degree to radian
R2D = 1/D2R;%                 convert radian to degree
Mug2mps2 = 9.80665E-6;%       convert micro-g to meter per second.^2
%% Run the Kalman filter + INS correction
%--------------------------------------------------------------------------
% Inputs:
%   rBEN_GPS        GPS estimated position (m)(@NED)
%   vBEN_GPS        GPS estimated velocity (m/s) (@NED)
%   tao_GPS         propagation interval (s)
%   L               previous latitude (rad)
%   lambda          previous longitude (rad)
%   h               previous height (m)
%   vBEN            prior estimated velocity @ NED frame (m/s)
%   TNB             prior estimated body to NED coordinate transformation matrix
%   IMU_bias_est    prior estimated IMU biases 
%   Cov_est         prior Kalman filter error covariance matrix
%   meas_fBIB       measured specific force of body w.r.t ECI  (m/s^2)
%   meas_wBIB       measured angular rate of body w.r.t ECI    (rad/s)
%   KF_unc          Kalman filter uncertainties
%   lGBB            Lever Arm
%--------------------------------------------------------------------------
% Outputs:
%   L_est           estimated latitude (rad)
%   lambda_est      estimated longitude (rad)
%   h_est           estimated height (m)
%   vBEN_est        estimated velocity @ NED frame (m/s)
%   TNB_est         estimated body to NED coordinate transformation matrix
%   IMU_bias_est    estimated IMU biases
%   Cov_est         estimated Kalman filter error covariance matrix
%==========================================================================
%% Preliminary calculations
%Constant Parameters
omega_ie = 7.292115E-5;   %Earth rotation rate (rad/s)
R0 = 6378137;        %WGS84 Equatorial radius (m)
e = 0.0818191908425; %WGS84 eccentricity

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details

% Begins

% Skew symmetric matrix of Earth rate
Omega_ie = Skew_symmetric([0,0,omega_ie]);
%--------------------------------------------------------------------------
%Calculate meridian & normal radius of the earth
[R_N,R_E]= Radii_of_curvature(L);
%abbreviation
REh = R_E + h;   RNh = R_N + h; 
vN  = est_v_eb_n(1);   vE = est_v_eb_n(2);   vD = est_v_eb_n(3);
sL  = sin(L);    cL = cos(L);    tL = tan(L);
%Calculate the earth rotation rate @ NED-frame
wNEN = [ vE / (REh);...
        -vN / (RNh);...
        -vE * tL / (REh)];
wEIN = omega_ie * [cL; 0; - sL];    
wNIN = wNEN + wEIN;
%Skew symmetric matrix
WEIN = Skew_symmetric(wEIN);
WNIN = Skew_symmetric(wNIN);
%==========================================================================
%% Determine state transition matrix using ??:
F11 = - WNIN;
F12 = [         0,   -1 / (REh),   0;...
        1 / (RNh),            0,   0;...
                0,   tL / (REh),   0];
%--------------------------------------------------------------------------                      
F13_31 = omega_ie * cL + vE / ((REh) * cL^2);                      
F13 = [ omega_ie * sL,   0,         vE / (REh)^2;...
               0,   0,        -vN / (RNh)^2;...
         F13_31 ,   0,   -vE * tL / (REh)^2];
%--------------------------------------------------------------------------
F21 = -Skew_symmetric(est_C_b_n_old * meas_fBIB);
%--------------------------------------------------------------------------
F22 = [ vD / RNh,    -2 * vE * tL / REh,   vN / RNh;...
               0,   (vN * tL + vD) / REh,  vE / REh;...
               0,                     0,          0];
F22 (2,1) = - F22 (1,2) /2; F22 (3,1) = - 2 * F22 (1,3); F22 (3,2) = - 2 * F22 (2,3);
F22 = F22 -2 * Skew_symmetric(wEIN);
%--------------------------------------------------------------------------
%geocentric radius using B.4:
rSEE = R0 / sqrt(1 - (e * sL)^2) * sqrt(cL^2 + (1 - e^2)^2 * sL^2);
%Calculate surface gravity using the Somigliana model using B.13, B.14:
sL2 = sL^2;
g0 = 9.7803253359 * (1 + 0.001931853 * sL2) / sqrt(1 - e^2 * sL2);
    %-------------------------------------------------------------
F23_1 = [                    - (vE / cL)^2 / REh - 2 * vE * omega_ie * cL;...
         vN * vE/ cL^2 / REh + 2 * vN * omega_ie * cL - 2 * vD * omega_ie * sL;... 
                                                   2 * vE * omega_ie * sL];
F23_3 = [         - (vE / REh)^2 * tL - vN * vD/ RNh^2;...
                     - (vN * vE * tL + vE * vD)/ REh^2;... 
          (vE / REh)^2 + (vN / RNh)^2 - 2 * g0 / rSEE];  
F23 = [ F23_1,   zeros(3,1),   F23_3];   
%--------------------------------------------------------------------------
F32 = [1 / RNh,               0,    0;...
             0,    1 / REh / cL,    0;...
             0,               0,   -1];
%--------------------------------------------------------------------------
F33 = [0, 0, -vN / RNh^2; vE * sL / REh / cL^2, 0, -vE / REh^2 / cL; 0, 0, 0];
%--------------------------------------------------------------------------           
% FN = zeros(15);
% FN(1:3,1:3)   = - WNIN ;
% FN(1:3,4:6)   = F12;
% FN(1:3,7:9)   = F13;
% FN(1:3,13:15) = est_C_b_n_old;
% FN(4:6,1:3)   = F21;
% FN(4:6,4:6)   = F22;
% FN(4:6,7:9)   = F23;
% FN(4:6,10:12) = est_C_b_n_old;
% FN(7:9,4:6)   = F32;
% FN(7:9,7:9)   = F33;
%--------------------------------------------------------------------------
%Continuous to Discrete Model: 
% delta_X'(t)  = FN * X(t)
% delta_X(k+1) = Phi_N * delta_X(k)
% [Phi_matrix,~] = Con2Dis(FN,0,tao_GPS,4);

% 1. Determine transition matrix using (14.50) (first-order approx)
Phi_matrix = eye(15);
Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) + F11 * tor_s;
Phi_matrix(1:3,4:6) = F12 * tor_s;;
Phi_matrix(1:3,7:9) = F13 * tor_s;;
Phi_matrix(1:3,13:15) = est_C_b_n_old * tor_s;

Phi_matrix(4:6,1:3) = F21 * tor_s;
Phi_matrix(4:6,4:6) = Phi_matrix(4:6,4:6) + F22 * tor_s;
Phi_matrix(4:6,7:9) = F23 * tor_s;
Phi_matrix(4:6,10:12) = est_C_b_n_old * tor_s;

Phi_matrix(7:9,4:6) = F32 * tor_s;
Phi_matrix(7:9,7:9) = Phi_matrix(7:9,7:9) + F33 * tor_s;

% Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) - Omega_ie * tor_s;
% Phi_matrix(1:3,13:15) = est_C_b_e_old * tor_s;
% Phi_matrix(4:6,1:3) = -tor_s * Skew_symmetric(est_C_b_e_old * meas_f_ib_b);
% Phi_matrix(4:6,4:6) = Phi_matrix(4:6,4:6) - 2 * Omega_ie * tor_s;
% geocentric_radius = R_0 / sqrt(1 - (e * sin(est_L_b_old))^2) *...
%     sqrt(cos(est_L_b_old)^2 + (1 - e^2)^2 * sin(est_L_b_old)^2); % from (2.137)
% Phi_matrix(4:6,7:9) = -tor_s * 2 * Gravity_ECEF(est_r_eb_e_old) /...
%     geocentric_radius * est_r_eb_e_old' / sqrt (est_r_eb_e_old' *...
%     est_r_eb_e_old);
% Phi_matrix(4:6,10:12) = est_C_b_e_old * tor_s;
% Phi_matrix(7:9,4:6) = eye(3) * tor_s;
%==========================================================================
%% Determine exact system noise covariance matrix using 1.44, 1.45:
tao_INS = 0.01;
    %init_ba_sd: Initial accel. bias uncertainty (m/s^2)
    init_ba_sd = 1000 * Mug2mps2;
    %accel_VRW: Accelerometer noise PSD (m^2 s^-3)
    accel_VRW = (1100 * Mug2mps2) ^ 2 * tao_INS;
    %accel_bias_PSD: Accelerometer bias random walk PSD (m^2 s^-5)
    accel_bias_PSD = 1.0E-5;
    %init_bg_sd: Initial gyro. bias uncertainty (rad/s)
    init_bg_sd = 20 * D2R / 3600;
    %gyro_ARW: Gyro noise PSD (rad^2/s)
    gyro_ARW = (1.1 * D2R / 60)^ 2 * tao_INS;
    %gyro_bias_PSD: Gyro bias random walk PSD (rad^2 s^-3)
    gyro_bias_PSD = 4.0E-11;
%-------------------------------------------------------------------------- 
%% IMU noise uncertainty
KF_unc.noise_IMU = [accel_VRW; gyro_ARW;  accel_bias_PSD;  gyro_bias_PSD;];

S_a  = KF_unc.noise_IMU (1);
S_g  = KF_unc.noise_IMU (2);
Sp_a = KF_unc.noise_IMU (3);
Sp_g = KF_unc.noise_IMU (4);
%Cartesian-to-curvilinear position change transformation matrix (p = T_prn * r_n)
T_prn = F32;
%--------------------------------------------------------------------------
Q11 = eye(3) * (S_g * tor_s + Sp_g * tor_s^3 / 3);
Q15 = Sp_g*tor_s^2*est_C_b_n_old / 2;
Q21 = (S_g*tor_s^2 / 2 + Sp_g*tor_s^4 / 4) * F21;
Q22 = (S_a*tor_s + Sp_a*tor_s^3/3)*eye(3) + (S_g*tor_s^3/3 + Sp_g*tor_s^5/5)*(F21*F21');
Q24 = Sp_a*tor_s^2*est_C_b_n_old / 2;
Q25 = Sp_g*tor_s^3*F21*est_C_b_n_old / 3;
Q31 = (S_g*tor_s^3 / 3 + Sp_g*tor_s^5 / 5) * T_prn * F21;
Q32 = (S_a*tor_s^2/2 + Sp_a*tor_s^4/4) * T_prn + (S_g*tor_s^4/4 + Sp_g*tor_s^6/6)*(T_prn*(F21*F21'));
Q33 = (S_a*tor_s^3/3 + Sp_a*tor_s^5/5)*T_prn*T_prn + (S_g*tor_s^5/5 + Sp_g*tor_s^7/7)*T_prn*(F21*F21')*T_prn;
Q34 = Sp_a*tor_s^3*T_prn*est_C_b_n_old/3;
Q35 = Sp_g*tor_s^4*T_prn*F21*est_C_b_n_old/4;
Q44 = Sp_a*tor_s*eye(3);
Q55 = Sp_g*tor_s*eye(3);
Q_prime_matrix = [Q11        Q21'    Q31'    zeros(3)    Q15;
       Q21        Q22     Q32'    Q24         Q25;
       Q31        Q32     Q33     Q34         Q35;
       zeros(3)   Q24'    Q34'    Q44         zeros(3);
       Q15'       Q25'    Q35'    zeros(3)    Q55;];
%==========================================================================
%% Propagation 
% 3. Propagate state estimates using (3.14) noting that all states are zero 
% due to closed-loop correction.
x_est_propagated(1:15,1) = 0;
%--------------------------------------------------------------------------
% 4. Propagate state estimation error covariance matrix using (3.46)
P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *...
    Phi_matrix' + 0.5 * Q_prime_matrix;
%==========================================================================
% MEASUREMENT UPDATE PHASE     

%rescale the latitude and longitude components
S_Ll = 1E3;
Sp = diag ([S_Ll, S_Ll, 1]);
% 5. Set-up measurement matrix using (14.115) 
H_matrix = zeros(6,15);
H_matrix(1:3,1:3)   =  Sp * (T_prn * Skew_symmetric(est_C_b_n_old * lGBB));
H_matrix(1:3,7:9)   = -Sp;
H_matrix(4:6,1:3)   = Skew_symmetric(est_C_b_n_old * cross(meas_wBIB,lGBB) - WEIN * est_C_b_n_old * lGBB);
H_matrix(4:6,4:6)   = -eye(3);
H_matrix(4:6,13:15) = est_C_b_n_old * Skew_symmetric(lGBB);

p = [L; lambda; h];
%Formulate measurement innovations using 1.53:
delta_z(1:3,1) = Sp * (GNSS_r_eb_n - p - T_prn * est_C_b_n_old * lGBB);
delta_z(4:6,1) = GNSS_v_eb_n - est_v_eb_n - est_C_b_n_old * cross(meas_wBIB,lGBB)+ WEIN * est_C_b_n_old * lGBB;
%--------------------------------------------------------------------------
%Set-up measurement matrix using 1.50, 1.52:
% H = zeros(6,15);
% H(1:3,1:3)   =  Sp * (T_prn * Skew_symmetric(est_C_b_n_old * lGBB));
% H(1:3,7:9)   = -Sp;
% H(4:6,1:3)   = Skew_symmetric(est_C_b_n_old * cross(meas_wBIB,lGBB) - WEIN * est_C_b_n_old * lGBB);
% H(4:6,4:6)   = -eye(3);
% H(4:6,13:15) = est_C_b_n_old * Skew_symmetric(lGBB);
%--------------------------------------------------------------------------
%Set-up measurement noise covariance matrix 
% R_matrix = KF_unc.Cov_meas;
% 6. Set-up measurement noise covariance matrix assuming all components of
% GNSS position and velocity are independent and have equal variance.
R_matrix(1:3,1:3) = eye(3) * LC_KF_config.pos_meas_SD^2;
R_matrix(1:3,4:6) = zeros(3);
R_matrix(4:6,1:3) = zeros(3);
R_matrix(4:6,4:6) = eye(3) * LC_KF_config.vel_meas_SD^2;
%==========================================================================
%% Estimation
% 7. Calculate Kalman gain using (3.21)
K_matrix = P_matrix_propagated * H_matrix' /(H_matrix *P_matrix_propagated * H_matrix' + R_matrix);
%--------------------------------------------------------------------------
% 9. Update state estimates using (3.24)
x_est_new = x_est_propagated + K_matrix * delta_z;
%--------------------------------------------------------------------------
% 10. Update state estimation error covariance matrix using (3.25)
P_matrix_new = (eye(15) - K_matrix * H_matrix) * P_matrix_propagated;
%==========================================================================
%% CLOSED-LOOP CORRECTION
% Correct attitude 
est_C_b_n_new = (eye(3) - Skew_symmetric(x_est_new(1:3))) * est_C_b_n_old;
%Orthonormalization
% est_C_b_n_new = OrthoNorm(est_C_b_n_new);
%--------------------------------------------------------------------------
%Correct velocity
est_v_eb_n_new = est_v_eb_n - x_est_new(4:6);
%--------------------------------------------------------------------------
%Correct position
L_est      = p(1) - x_est_new(7);
lambda_est = p(2) - x_est_new(8);
h_est      = p(3) - x_est_new(9);
%--------------------------------------------------------------------------
% Update IMU bias estimates
est_IMU_bias_new = est_IMU_bias_old + x_est_new(10:15);
%==========================================================================
end