%8/8/17 Script developed to curate the data into a form more usable by
%error_11h
%8/16/18 Script modified for LengthScaleFactor (LSF)
%10/17/18 Script modified for size extension (LE)
%1/11/19 Script modified for Model Predictive Control (MPC)
    %1/16/19 Script corrected previous indexing error.
    
%v11b is for horizontal aggressive maneuver
%v11c is for aggressive maneuver
%v12h is for sum of prime number sines

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% Assign moth variables IF NECESSARY
% disp('Assign moth parameters')
L1 = 0.908;
L2 = 1.7475;

bl = 5.311; %The body length of the organism in cm - THIS WILL CHANGE
          %FOR DAUBER (2.31), BEE (1.153), MOTH (5.311).

LSF_1 = 1;
          
%% List the directory stuff
listOfconLSF_1_LE_0_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_0.mat'); 
listOfconLSF_1_LE_p2_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p2.mat'); 
listOfconLSF_1_LE_p4_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p4.mat'); 
listOfconLSF_1_LE_p6_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p6.mat'); 
listOfconLSF_1_LE_p8_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p8.mat'); 
listOfconLSF_1_LE_1_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_1.mat'); 
listOfconLSF_1_LE_2_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_2.mat'); 
%the asterisk is a wildcard
%The dir function returns a "listing" of an M x 1 "structure." The 
%structure has five fields in this case listing: name, date, byte, isdir, 
%datenum.
%I used the wildcard because I know the number of text files will
%definitely increase as we gather more data.
%For more information enter   help dir   into MATLAB mainframe
disp('Load time vector')
load('../SimData_MPC/Tstore_MPC_hws_sp.mat')

%% Assign the imported numbers to internal variables
disp('Assign internal variables')

hws = 500; 
timestep = numel(Tstore(1:(end-100)))/hws;
tderiv = Tstore(2);
t_end = Tstore((end-100+1));

%Create the sum of primes signal
signal_amp = 5; %in cm 
prime_f = [0.2, 0.3, 0.5, 0.7, 1.1, 1.7, 2.9, 4.3, 7.9, 13.7, 19.9]; %in Hz
prime_a = (signal_amp./(2*pi.*prime_f)).*(2*pi.*prime_f(1)); %in cm
prime_ph = zeros(1,numel(prime_f)); %Not sure if we'll require a phase

%Goal criteria (for the cost function)
y_g = prime_a(1)*sin(2*pi*prime_f(1)*Tstore(1:(end-100+1)) + prime_ph(1)) +...
    prime_a(2)*sin(2*pi*prime_f(2)*Tstore(1:(end-100+1)) + prime_ph(2)) +...
    prime_a(3)*sin(2*pi*prime_f(3)*Tstore(1:(end-100+1)) + prime_ph(3)) +...
    prime_a(4)*sin(2*pi*prime_f(4)*Tstore(1:(end-100+1)) + prime_ph(4)) +...
    prime_a(5)*sin(2*pi*prime_f(5)*Tstore(1:(end-100+1)) + prime_ph(5)) +...
    prime_a(6)*sin(2*pi*prime_f(6)*Tstore(1:(end-100+1)) + prime_ph(6)) +...
    prime_a(7)*sin(2*pi*prime_f(7)*Tstore(1:(end-100+1)) + prime_ph(7)) +...
    prime_a(8)*sin(2*pi*prime_f(8)*Tstore(1:(end-100+1)) + prime_ph(8)) +...
    prime_a(9)*sin(2*pi*prime_f(9)*Tstore(1:(end-100+1)) + prime_ph(9)) +...
    prime_a(10)*sin(2*pi*prime_f(10)*Tstore(1:(end-100+1)) + prime_ph(10)) +...
    prime_a(11)*sin(2*pi*prime_f(11)*Tstore(1:(end-100+1)) + prime_ph(11));

ydot_g = 2*pi*prime_a(1)*prime_f(1)*cos(2*pi*prime_f(1)*Tstore(1:(end-100+1)) + prime_ph(1)) +...
    2*pi*prime_a(2)*prime_f(2)*cos(2*pi*prime_f(2)*Tstore(1:(end-100+1)) + prime_ph(2)) +...
    2*pi*prime_a(3)*prime_f(3)*cos(2*pi*prime_f(3)*Tstore(1:(end-100+1)) + prime_ph(3)) +...
    2*pi*prime_a(4)*prime_f(4)*cos(2*pi*prime_f(4)*Tstore(1:(end-100+1)) + prime_ph(4)) +...
    2*pi*prime_a(5)*prime_f(5)*cos(2*pi*prime_f(5)*Tstore(1:(end-100+1)) + prime_ph(5)) +...
    2*pi*prime_a(6)*prime_f(6)*cos(2*pi*prime_f(6)*Tstore(1:(end-100+1)) + prime_ph(6)) +...
    2*pi*prime_a(7)*prime_f(7)*cos(2*pi*prime_f(7)*Tstore(1:(end-100+1)) + prime_ph(7)) +...
    2*pi*prime_a(8)*prime_f(8)*cos(2*pi*prime_f(8)*Tstore(1:(end-100+1)) + prime_ph(8)) +...
    2*pi*prime_a(9)*prime_f(9)*cos(2*pi*prime_f(9)*Tstore(1:(end-100+1)) + prime_ph(9)) +...
    2*pi*prime_a(10)*prime_f(10)*cos(2*pi*prime_f(10)*Tstore(1:(end-100+1)) + prime_ph(10)) +...
    2*pi*prime_a(11)*prime_f(11)*cos(2*pi*prime_f(11)*Tstore(1:(end-100+1)) + prime_ph(11));

y_g = y_g';
x_g = zeros(1,numel(y_g));
theta_g = zeros(1,numel(y_g));
theta_g(1,:) = pi/4;

ydot_g = ydot_g';
xdot_g = zeros(1,numel(y_g));
thetadot_g = zeros(1,numel(y_g));

NumOfconLSF_1_LE_0_Files = numel(listOfconLSF_1_LE_0_Files);
NumOfconLSF_1_LE_p2_Files = numel(listOfconLSF_1_LE_p2_Files);
NumOfconLSF_1_LE_p4_Files = numel(listOfconLSF_1_LE_p4_Files);
NumOfconLSF_1_LE_p6_Files = numel(listOfconLSF_1_LE_p6_Files);
NumOfconLSF_1_LE_p8_Files = numel(listOfconLSF_1_LE_p8_Files);
NumOfconLSF_1_LE_1_Files = numel(listOfconLSF_1_LE_1_Files);
NumOfconLSF_1_LE_2_Files = numel(listOfconLSF_1_LE_2_Files);

%% Import all the variables
disp('Importing curated data')

%con - LSF_1_LE_0
load('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_0.mat')

%con - LSF_1_LE_p2
load('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p2.mat')

%con - LSF_1_LE_p4
load('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p4.mat')

%con - LSF_1_LE_p6
load('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p6.mat')

%con - LSF_1_LE_p8
load('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p8.mat')

%con - LSF_1_LE_1
load('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_1.mat')

%con - LSF_1_LE_2
load('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_2.mat')
load('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_2.mat')

%numel
load('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_2.mat')

%imported cost
load('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_2.mat')

load('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_2.mat')

%tepfr -- tracking error per full run
load('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_2.mat')

%Mean of cost per full run -- for stats purposes
load('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_2.mat')

%Tracking error w.r.t. time
load('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_2.mat')

disp('Curated data is loaded')

%% Find the min & max values for...
%x
xMax_vec = [max(max(mean_x_con_LSF_1_LE_0 + std_x_con_LSF_1_LE_0)),...
    max(max(mean_x_con_LSF_1_LE_p2 + std_x_con_LSF_1_LE_p2)),...
    max(max(mean_x_con_LSF_1_LE_p4 + std_x_con_LSF_1_LE_p4)),...
    max(max(mean_x_con_LSF_1_LE_p6 + std_x_con_LSF_1_LE_p6)),...
    max(max(mean_x_con_LSF_1_LE_p8 + std_x_con_LSF_1_LE_p8)),...
    max(max(mean_x_con_LSF_1_LE_1 + std_x_con_LSF_1_LE_1)),...
    max(max(mean_x_con_LSF_1_LE_2 + std_x_con_LSF_1_LE_2))];
x_max = max(xMax_vec);

xMin_vec = [min(min(mean_x_con_LSF_1_LE_0 - std_x_con_LSF_1_LE_0)),...
    min(min(mean_x_con_LSF_1_LE_p2 - std_x_con_LSF_1_LE_p2)),...
    min(min(mean_x_con_LSF_1_LE_p4 - std_x_con_LSF_1_LE_p4)),...
    min(min(mean_x_con_LSF_1_LE_p6 - std_x_con_LSF_1_LE_p6)),...
    min(min(mean_x_con_LSF_1_LE_p8 - std_x_con_LSF_1_LE_p8)),...
    min(min(mean_x_con_LSF_1_LE_1 - std_x_con_LSF_1_LE_1)),...
    min(min(mean_x_con_LSF_1_LE_2 - std_x_con_LSF_1_LE_2))];
x_min = min(xMin_vec);

%y
yMax_vec = [max(max(mean_y_con_LSF_1_LE_0 + std_y_con_LSF_1_LE_0)),...
    max(max(mean_y_con_LSF_1_LE_p2 + std_y_con_LSF_1_LE_p2)),...
    max(max(mean_y_con_LSF_1_LE_p4 + std_y_con_LSF_1_LE_p4)),...
    max(max(mean_y_con_LSF_1_LE_p6 + std_y_con_LSF_1_LE_p6)),...
    max(max(mean_y_con_LSF_1_LE_p8 + std_y_con_LSF_1_LE_p8)),...
    max(max(mean_y_con_LSF_1_LE_1 + std_y_con_LSF_1_LE_1)),...
    max(max(mean_y_con_LSF_1_LE_2 + std_y_con_LSF_1_LE_2))];
y_max = max(yMax_vec);

yMin_vec = [min(min(mean_y_con_LSF_1_LE_0 - std_y_con_LSF_1_LE_0)),...
    min(min(mean_y_con_LSF_1_LE_p2 - std_y_con_LSF_1_LE_p2)),...
    min(min(mean_y_con_LSF_1_LE_p4 - std_y_con_LSF_1_LE_p4)),...
    min(min(mean_y_con_LSF_1_LE_p6 - std_y_con_LSF_1_LE_p6)),...
    min(min(mean_y_con_LSF_1_LE_p8 - std_y_con_LSF_1_LE_p8)),...
    min(min(mean_y_con_LSF_1_LE_1 - std_y_con_LSF_1_LE_1)),...
    min(min(mean_y_con_LSF_1_LE_2 - std_y_con_LSF_1_LE_2))];
y_min = min(yMin_vec);

%theta
thetaMax_vec = [max(max(mean_theta_con_LSF_1_LE_0 + std_theta_con_LSF_1_LE_0)),...
    max(max(mean_theta_con_LSF_1_LE_p2 + std_theta_con_LSF_1_LE_p2)),...
    max(max(mean_theta_con_LSF_1_LE_p4 + std_theta_con_LSF_1_LE_p4)),...
    max(max(mean_theta_con_LSF_1_LE_p6 + std_theta_con_LSF_1_LE_p6)),...
    max(max(mean_theta_con_LSF_1_LE_p8 + std_theta_con_LSF_1_LE_p8)),...
    max(max(mean_theta_con_LSF_1_LE_1 + std_theta_con_LSF_1_LE_1)),...
    max(max(mean_theta_con_LSF_1_LE_2 + std_theta_con_LSF_1_LE_2))];
theta_max = max(thetaMax_vec);

thetaMin_vec = [min(min(mean_theta_con_LSF_1_LE_0 - std_theta_con_LSF_1_LE_0)),...
    min(min(mean_theta_con_LSF_1_LE_p2 - std_theta_con_LSF_1_LE_p2)),...
    min(min(mean_theta_con_LSF_1_LE_p4 - std_theta_con_LSF_1_LE_p4)),...
    min(min(mean_theta_con_LSF_1_LE_p6 - std_theta_con_LSF_1_LE_p6)),...
    min(min(mean_theta_con_LSF_1_LE_p8 - std_theta_con_LSF_1_LE_p8)),...
    min(min(mean_theta_con_LSF_1_LE_1 - std_theta_con_LSF_1_LE_1)),...
    min(min(mean_theta_con_LSF_1_LE_2 - std_theta_con_LSF_1_LE_2))];
theta_min = min(thetaMin_vec);

%phi
phiMax_vec = [max(max(mean_phi_con_LSF_1_LE_0 + std_phi_con_LSF_1_LE_0)),...
    max(max(mean_phi_con_LSF_1_LE_p2 + std_phi_con_LSF_1_LE_p2)),...
    max(max(mean_phi_con_LSF_1_LE_p4 + std_phi_con_LSF_1_LE_p4)),...
    max(max(mean_phi_con_LSF_1_LE_p6 + std_phi_con_LSF_1_LE_p6)),...
    max(max(mean_phi_con_LSF_1_LE_p8 + std_phi_con_LSF_1_LE_p8)),...
    max(max(mean_phi_con_LSF_1_LE_1 + std_phi_con_LSF_1_LE_1)),...
    max(max(mean_phi_con_LSF_1_LE_2 + std_phi_con_LSF_1_LE_2))];
phi_max = max(phiMax_vec);

phiMin_vec = [min(min(mean_phi_con_LSF_1_LE_0 - std_phi_con_LSF_1_LE_0)),...
    min(min(mean_phi_con_LSF_1_LE_p2 - std_phi_con_LSF_1_LE_p2)),...
    min(min(mean_phi_con_LSF_1_LE_p4 - std_phi_con_LSF_1_LE_p4)),...
    min(min(mean_phi_con_LSF_1_LE_p6 - std_phi_con_LSF_1_LE_p6)),...
    min(min(mean_phi_con_LSF_1_LE_p8 - std_phi_con_LSF_1_LE_p8)),...
    min(min(mean_phi_con_LSF_1_LE_1 - std_phi_con_LSF_1_LE_1)),...
    min(min(mean_phi_con_LSF_1_LE_2 - std_phi_con_LSF_1_LE_2))];
phi_min = min(phiMin_vec);

%beta
betaMax_vec = [max(max(mean_beta_con_LSF_1_LE_0 + std_beta_con_LSF_1_LE_0)),...
    max(max(mean_beta_con_LSF_1_LE_p2 + std_beta_con_LSF_1_LE_p2)),...
    max(max(mean_beta_con_LSF_1_LE_p4 + std_beta_con_LSF_1_LE_p4)),...
    max(max(mean_beta_con_LSF_1_LE_p6 + std_beta_con_LSF_1_LE_p6)),...
    max(max(mean_beta_con_LSF_1_LE_p8 + std_beta_con_LSF_1_LE_p8)),...
    max(max(mean_beta_con_LSF_1_LE_1 + std_beta_con_LSF_1_LE_1)),...
    max(max(mean_beta_con_LSF_1_LE_2 + std_beta_con_LSF_1_LE_2))];
beta_max = max(betaMax_vec);

betaMin_vec = [min(min(mean_beta_con_LSF_1_LE_0 - std_beta_con_LSF_1_LE_0)),...
    min(min(mean_beta_con_LSF_1_LE_p2 - std_beta_con_LSF_1_LE_p2)),...
    min(min(mean_beta_con_LSF_1_LE_p4 - std_beta_con_LSF_1_LE_p4)),...
    min(min(mean_beta_con_LSF_1_LE_p6 - std_beta_con_LSF_1_LE_p6)),...
    min(min(mean_beta_con_LSF_1_LE_p8 - std_beta_con_LSF_1_LE_p8)),...
    min(min(mean_beta_con_LSF_1_LE_1 - std_beta_con_LSF_1_LE_1)),...
    min(min(mean_beta_con_LSF_1_LE_2 - std_beta_con_LSF_1_LE_2))];
beta_min = min(betaMin_vec);

%tepfr
tepfrMax_vec = [max(max(mean_tepfr_con_LSF_1_LE_0)),...
    max(max(mean_tepfr_con_LSF_1_LE_p2)),...
    max(max(mean_tepfr_con_LSF_1_LE_p4)),...
    max(max(mean_tepfr_con_LSF_1_LE_p6)),...
    max(max(mean_tepfr_con_LSF_1_LE_p8)),...
    max(max(mean_tepfr_con_LSF_1_LE_1)),...
    max(max(mean_tepfr_con_LSF_1_LE_2))];
tepfr_max = max(tepfrMax_vec);

tepfrMin_vec = [min(min(mean_tepfr_con_LSF_1_LE_0)),...
    min(min(mean_tepfr_con_LSF_1_LE_p2)),...
    min(min(mean_tepfr_con_LSF_1_LE_p4)),...
    min(min(mean_tepfr_con_LSF_1_LE_p6)),...
    min(min(mean_tepfr_con_LSF_1_LE_p8)),...
    min(min(mean_tepfr_con_LSF_1_LE_1)),...
    min(min(mean_tepfr_con_LSF_1_LE_2))];
tepfr_min = min(tepfrMin_vec);

%cost
costMax_vec = [max(max(mean_cost_con_LSF_1_LE_0)),...
    max(max(mean_cost_con_LSF_1_LE_p2)),...
    max(max(mean_cost_con_LSF_1_LE_p4)),...
    max(max(mean_cost_con_LSF_1_LE_p6)),...
    max(max(mean_cost_con_LSF_1_LE_p8)),...
    max(max(mean_cost_con_LSF_1_LE_1)),...
    max(max(mean_cost_con_LSF_1_LE_2))];
cost_max = max(costMax_vec);

costMin_vec = [min(min(mean_cost_con_LSF_1_LE_0)),...
    min(min(mean_cost_con_LSF_1_LE_p2)),...
    min(min(mean_cost_con_LSF_1_LE_p4)),...
    min(min(mean_cost_con_LSF_1_LE_p6)),...
    min(min(mean_cost_con_LSF_1_LE_p8)),...
    min(min(mean_cost_con_LSF_1_LE_1)),...
    min(min(mean_cost_con_LSF_1_LE_2))];
cost_min = min(costMin_vec);

%% Set up internal values for Figure 1

abdoflip(1:numel(Tstore)) = 2*pi;
livebetarange(1:numel(Tstore)) = 20*(pi/180); 

%% Control variables w.r.t. time
fig1 = figure(1);

%x - con LSF_1_LE_0
subplot(5,7,1)
% subplot(5,19,10)
Fig1_1 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_con_LSF_1_LE_0,std_x_con_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_1.mainLine.LineWidth = 2;
title({['LSF: ',num2str(LSF_1)];'LE: 0'})
ylabel({'x';'(cm)'});
% axis([0, Tstore(end-100+1), x_min, x_max])

%x - con LSF_1_LE_p2
subplot(5,7,2)
Fig1_2 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_con_LSF_1_LE_p2,...
    std_x_con_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_2.mainLine.LineWidth = 2;
title({['LSF: ',num2str(LSF_1)];'LE: 0.2'})
% axis([0, Tstore(end-100+1), x_min, x_max])

%x - con LSF_1_LE_p4
subplot(5,7,3)
Fig1_3 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_con_LSF_1_LE_p4,...
    std_x_con_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_3.mainLine.LineWidth = 2;
title({['LSF: ',num2str(LSF_1)];'LE: 0.4'})
% axis([0, Tstore(end-100+1), x_min, x_max])

%x - con LSF_1_LE_p6
subplot(5,7,4)
Fig1_4 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_con_LSF_1_LE_p6,...
    std_x_con_LSF_1_LE_p6,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_4.mainLine.LineWidth = 2;
title({['LSF: ',num2str(LSF_1)];'LE: 0.6'})
% axis([0, Tstore(end-100+1), x_min, x_max])

%x - con LSF_1_LE_p8
subplot(5,7,5)
Fig1_5 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_con_LSF_1_LE_p8,...
    std_x_con_LSF_1_LE_p8,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_5.mainLine.LineWidth = 2;
title({['LSF: ',num2str(LSF_1)];'LE: 0.8'})
% axis([0, Tstore(end-100+1), x_min, x_max])

%x - con LSF_1_LE_1
subplot(5,7,6)
Fig1_6 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_con_LSF_1_LE_1,...
    std_x_con_LSF_1_LE_1,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_6.mainLine.LineWidth = 2;
title({['LSF: ',num2str(LSF_1)];'LE: 1'})
% axis([0, Tstore(end-100+1), x_min, x_max])

%x - con LSF_1_LE_2
subplot(5,7,7)
Fig1_7 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_con_LSF_1_LE_2,...
    std_x_con_LSF_1_LE_2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_7.mainLine.LineWidth = 2;
title({['LSF: ',num2str(LSF_1)];'LE: 2'})
% axis([0, Tstore(end-100+1), x_min, x_max])

%y - con LSF_1_LE_0
subplot(5,7,8)
Fig1_8 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_con_LSF_1_LE_0,...
    std_y_con_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_8.mainLine.LineWidth = 1.2; 
hold on;
ylabel({'y';'(cm)'});
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])

%y - con LSF_1_LE_p2
subplot(5,7,9)
Fig1_9 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_con_LSF_1_LE_p2,...
    std_y_con_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_9.mainLine.LineWidth = 1.2; 
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])

%y - con LSF_1_LE_p4
subplot(5,7,10)
Fig1_10 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_con_LSF_1_LE_p4,...
    std_y_con_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_10.mainLine.LineWidth = 1.2;
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])

%y - con LSF_1_LE_p6
subplot(5,7,11)
Fig1_11 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_con_LSF_1_LE_p6,...
    std_y_con_LSF_1_LE_p6,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_11.mainLine.LineWidth = 1.2; 
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])

%y - con LSF_1_LE_p8
subplot(5,7,12)
Fig1_12 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_con_LSF_1_LE_p8,...
    std_y_con_LSF_1_LE_p8,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_12.mainLine.LineWidth = 1.2;
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])

%y - con LSF_1_LE_1
subplot(5,7,13)
Fig1_13 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_con_LSF_1_LE_1,...
    std_y_con_LSF_1_LE_1,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_13.mainLine.LineWidth = 1.2;
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])

%y - con LSF_1_LE_2
subplot(5,7,14)
Fig1_14 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_con_LSF_1_LE_2,...
    std_y_con_LSF_1_LE_2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_14.mainLine.LineWidth = 1.2;
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])

%theta - con LSF_1_LE_0
subplot(5,7,15)
Fig1_15 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_con_LSF_1_LE_0,...
    (180/pi)*std_theta_con_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
ylabel({'theta';'(degrees)'});
Fig1_15.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])

%theta - con LSF_1_LE_p2
subplot(5,7,16)
Fig1_16 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_con_LSF_1_LE_p2,...
    (180/pi)*std_theta_con_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_16.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])

%theta - con LSF_1_LE_p4
subplot(5,7,17)
Fig1_17 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_con_LSF_1_LE_p4,...
    (180/pi)*std_theta_con_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_17.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])

%theta - con LSF_1_LE_p6
subplot(5,7,18)
Fig1_18 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_con_LSF_1_LE_p6,...
    (180/pi)*std_theta_con_LSF_1_LE_p6,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_18.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])

%theta - con LSF_1_LE_p8
subplot(5,7,19)
Fig1_19 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_con_LSF_1_LE_p8,...
    (180/pi)*std_theta_con_LSF_1_LE_p8,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_19.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])

%theta - con LSF_1_LE_1
subplot(5,7,20)
Fig1_20 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_con_LSF_1_LE_1,...
    (180/pi)*std_theta_con_LSF_1_LE_1,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_20.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])

%theta - con LSF_1_LE_2
subplot(5,7,21)
Fig1_21 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_con_LSF_1_LE_2,...
    (180/pi)*std_theta_con_LSF_1_LE_2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_21.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])

%phi - con LSF_1_LE_0
subplot(5,7,22)
Fig1_22 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_con_LSF_1_LE_0,...
    (180/pi)*std_phi_con_LSF_1_LE_0,...
    'lineprops', '-k','transparent',false,'patchSaturation',0.075);
Fig1_22.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
ylabel({'phi';'(degrees)'});
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])

%phi - con LSF_1_LE_p2
subplot(5,7,23)
Fig1_23 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_con_LSF_1_LE_p2,...
    (180/pi)*std_phi_con_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_23.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])

%phi - con LSF_1_LE_p4
subplot(5,7,24)
Fig1_24 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_con_LSF_1_LE_p4,...
    (180/pi)*std_phi_con_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_24.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])

%phi - con LSF_1_LE_p6
subplot(5,7,25)
Fig1_25 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_con_LSF_1_LE_p6,...
    (180/pi)*std_phi_con_LSF_1_LE_p6,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig1_25.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])

%phi - con LSF_1_LE_p8
subplot(5,7,26)
Fig1_26 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_con_LSF_1_LE_p8,...
    (180/pi)*std_phi_con_LSF_1_LE_p8,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_26.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])

%phi - con LSF_1_LE_1
subplot(5,7,27)
Fig1_27 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_con_LSF_1_LE_1,...
    (180/pi)*std_phi_con_LSF_1_LE_1,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig1_27.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])

%phi - con LSF_1_LE_2
subplot(5,7,28)
Fig1_28 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_con_LSF_1_LE_2,...
    (180/pi)*std_phi_con_LSF_1_LE_2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig1_28.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])

%beta - con LSF_1_LE_0
subplot(5,7,29)
Fig1_29 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_con_LSF_1_LE_0,...
    (180/pi)*std_beta_con_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_29.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfconLSF_1_LE_0_Files),')']}); 
ylabel({'beta';'(degrees)'});
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);

%beta - con LSF_1_LE_p2
subplot(5,7,30)
Fig1_30 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_con_LSF_1_LE_p2,...
    (180/pi)*std_beta_con_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_30.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);

%beta - con LSF_1_LE_p4
subplot(5,7,31)
Fig1_31 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_con_LSF_1_LE_p4,...
    (180/pi)*std_beta_con_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig1_31.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);

%beta - con LSF_1_LE_p6
subplot(5,7,32)
Fig1_32 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_con_LSF_1_LE_p6,...
    (180/pi)*std_beta_con_LSF_1_LE_p6,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig1_32.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',-120:40:180);

%beta - con LSF_1_LE_p8
subplot(5,7,33)
Fig1_33 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_con_LSF_1_LE_p8,...
    (180/pi)*std_beta_con_LSF_1_LE_p8,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig1_33.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',-120:40:180);

%beta - con LSF_1_LE_1
subplot(5,7,34)
Fig1_34 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_con_LSF_1_LE_1,...
    (180/pi)*std_beta_con_LSF_1_LE_1,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig1_34.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfconLSF_1_LE_1_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',-120:40:180);

%beta - con LSF_1_LE_2
subplot(5,7,35)
Fig1_35 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_con_LSF_1_LE_2,...
    (180/pi)*std_beta_con_LSF_1_LE_2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig1_35.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfconLSF_1_LE_2_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',-120:40:180);

set(fig1, 'Position', [0,0, 1250, 2500])

saveas(gcf,'../Figures_MPC/LSF_1/Fig1_MPC_statevars_time_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig1_MPC_statevars_time_LSF_1_v12h.jpg');

disp('Figure 1 done')

%% Tracking error w.r.t. time

figure(2);
%y motion- goal
subplot(2,1,1)
% subplot(8,1,1)
plot(Tstore(1:(end-100+1)),y_g,'-k','LineWidth',2);
ylabel({'Goal y-motion'; '(cm)'}); 
% xlabel({'0'; ['(n = ',num2str(NumOfconLSF_1_LE_0_Files),')']});
axis([0, Tstore(end-100+1), -inf,inf])

%Tracking error - con LSF_1_LE_0
subplot(2,1,2)
% subplot(8,1,2)
Fig2_1 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_con_LSF_1_LE_0,...
    std_dist_con_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_1.mainLine.LineWidth = 2;
ylabel({'Tracking error'; '(cm)'}); 
xlabel({'LE: 0'; ['(n = ',num2str(NumOfconLSF_1_LE_0_Files),')']});
axis([0, Tstore(end-100+1), -inf,inf])
% set(gca,'ytick',0:0.1:0.6);

% %Tracking error - con LSF_1_LE_p2
% subplot(8,1,3)
% Fig2_2 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_con_LSF_1_LE_p2,...
%     std_dist_con_LSF_1_LE_p2,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_2.mainLine.LineWidth = 2;
% ylabel({'Tracking error'; '(cm)'}); 
% xlabel({'0.2'; ['(n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']});
% axis([0, Tstore(end-100+1), -inf,inf])
% 
% %Tracking error - con LSF_1_LE_p4
% subplot(8,1,4)
% Fig2_3 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_con_LSF_1_LE_p4,...
%     std_dist_con_LSF_1_LE_p4,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_3.mainLine.LineWidth = 2;
% ylabel({'Tracking error'; '(cm)'}); 
% xlabel({'0.4'; ['(n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']});
% axis([0, Tstore(end-100+1), -inf,inf])
% 
% %Tracking error - con LSF_1_LE_p6
% subplot(8,1,5)
% Fig2_4 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_con_LSF_1_LE_p6,...
%     std_dist_con_LSF_1_LE_p6,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_4.mainLine.LineWidth = 2;
% ylabel({'Tracking error'; '(cm)'}); 
% xlabel({'0.6'; ['(n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']});
% axis([0, Tstore(end-100+1), -inf,inf])
% 
% %Tracking error - con LSF_1_LE_p8
% subplot(8,1,6)
% Fig2_5 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_con_LSF_1_LE_p8,...
%     std_dist_con_LSF_1_LE_p8,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_5.mainLine.LineWidth = 2;
% ylabel({'Tracking error'; '(cm)'}); 
% xlabel({'0.8'; ['(n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']});
% axis([0, Tstore(end-100+1), -inf,inf])

% %Tracking error - con LSF_1_LE_1
% subplot(8,1,7)
% Fig2_6 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_con_LSF_1_LE_1,...
%     std_dist_con_LSF_1_LE_1,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_6.mainLine.LineWidth = 2;
% ylabel({'Tracking error'; '(cm)'}); 
% xlabel({'1'; ['(n = ',num2str(NumOfconLSF_1_LE_1_Files),')']});
% axis([0, Tstore(end-100+1), -inf,inf])
% 
% %Tracking error - con LSF_1_LE_2
% subplot(8,1,8)
% Fig2_7 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_con_LSF_1_LE_2,...
%     std_dist_con_LSF_1_LE_2,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_7.mainLine.LineWidth = 2;
% ylabel({'Tracking error'; '(cm)'}); 
% xlabel({'2'; ['(n = ',num2str(NumOfconLSF_1_LE_2_Files),')',' Time (sec)']});
% axis([0, Tstore(end-100+1), -inf,inf])

saveas(gcf,'../Figures_MPC/LSF_1/Fig2_MPC_error_time_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig2_MPC_error_time_LSF_1_v12h.jpg');

disp('Figure 2 done')

%% Tracking error violin plots unscaled

mothL = 2*L1 + 2*L2; %in cm

fig3 = figure(3);
%Tracking error - con LSF_1_LE_0
subplot(1,7,1)
distributionPlot(mean_tepfr_con_LSF_1_LE_0,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_0_Files),')'])
ylabel('Mean of tracking error (cm)')
% axis([0, 2, tepfr_min, tepfr_max])
% set(gca,'ytick',0:0.5:1);

%Tracking error - con LSF_1_LE_p2
subplot(1,7,2) 
distributionPlot(mean_tepfr_con_LSF_1_LE_p2,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.2'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p2_Files),')'])
% axis([0, 2, tepfr_min, tepfr_max])
% set(gca,'ytick',0:0.5:1);

%Tracking error - con LSF_1_LE_p4
subplot(1,7,3)
distributionPlot(mean_tepfr_con_LSF_1_LE_p4,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.4'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p4_Files),')'])
% axis([0, 2, tepfr_min, tepfr_max])
% set(gca,'ytick',0:0.5:1);

%Tracking error - con LSF_1_LE_p6
subplot(1,7,4)
distributionPlot(mean_tepfr_con_LSF_1_LE_p6,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.6'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p6_Files),')'])
% axis([0, 2, tepfr_min, tepfr_max])
% set(gca,'ytick',0:10:80);

%Tracking error - con LSF_1_LE_p8
subplot(1,7,5)
distributionPlot(mean_tepfr_con_LSF_1_LE_p8,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.8'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p8_Files),')'])
% axis([0, 2, tepfr_min, tepfr_max])
% set(gca,'ytick',0:10:80);

%Tracking error - con LSF_1_LE_1
subplot(1,7,6)
distributionPlot(mean_tepfr_con_LSF_1_LE_1,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 1'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_1_Files),')'])
% axis([0, 2, tepfr_min, tepfr_max])
% set(gca,'ytick',0:10:80);

%Tracking error - con LSF_1_LE_2
subplot(1,7,7)
distributionPlot(mean_tepfr_con_LSF_1_LE_2,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 2'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_2_Files),')'])
% axis([0, 2, tepfr_min, tepfr_max])
% set(gca,'ytick',0:10:80);

set(fig3, 'Position', [120,120,850,525])

saveas(gcf,'../Figures_MPC/LSF_1/Fig3_MPC_TrackErrorViolin_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig3_MPC_TrackErrorVoilin_LSF_1_v12h.jpg');

disp('Figure 3 done') 

%% Tracking error scaled

fig4 = figure(4);
%Tracking error - con LSF_1_LE_0
subplot(1,7,1)
distributionPlot(mean_tepfr_con_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_0_Files),')'])
ylabel({'Mean of tracking error';'normalized to body length'})
% axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
% set(gca,'ytick',0:0.05:0.2);

%Tracking error - con LSF_1_LE_p2
subplot(1,7,2)
distributionPlot(mean_tepfr_con_LSF_1_LE_p2./(mothL*LSF_1),'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.2'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p2_Files),')'])
% axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
% set(gca,'ytick',0:0.05:0.2);

%Tracking error - con LSF_1_LE_p4
subplot(1,7,3)
distributionPlot(mean_tepfr_con_LSF_1_LE_p4./(mothL*LSF_1),'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.4'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p4_Files),')'])
% axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
% set(gca,'ytick',0:0.05:0.2);

%Tracking error - con LSF_1_LE_p6
subplot(1,7,4)
distributionPlot(mean_tepfr_con_LSF_1_LE_p6./(mothL*LSF_1),'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.6'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p6_Files),')'])
% axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
% set(gca,'ytick',0:2:16);

%Tracking error - con LSF_1_LE_p8
subplot(1,7,5)
distributionPlot(mean_tepfr_con_LSF_1_LE_p8./(mothL*LSF_1),'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.8'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p8_Files),')'])
% axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
% set(gca,'ytick',0:2:16);

%Tracking error - con LSF_1_LE_1
subplot(1,7,6)
distributionPlot(mean_tepfr_con_LSF_1_LE_1./(mothL*LSF_1),'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 1'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_1_Files),')'])
% axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
% set(gca,'ytick',0:2:16);

%Tracking error - con LSF_1_LE_2
subplot(1,7,7)
distributionPlot(mean_tepfr_con_LSF_1_LE_2./(mothL*LSF_1),'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 2'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_2_Files),')'])
% axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
% set(gca,'ytick',0:2:16);

set(fig4, 'Position', [120,120,850,525])

saveas(gcf,'../Figures_MPC/LSF_1/Fig4_MPC_TrackErrorViolin_lengthscaled_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig4_MPC_TrackErrorVoilin_lengthscaled_LSF_1_v12h.jpg');

disp('Figure 4 done') 

%% Number of elements w.r.t. time

% figure(5);
% %Number of elements - con LSF_1_LE_0
% subplot(7,1,1)
% plot(Tstore, numel_con_LSF_1_LE_0,'k-','LineWidth',2)
% xlabel(['Time (s), con - LSF: 1, (n = ',...
%     num2str(NumOfconLSF_1_LE_0_Files),') at start, ',...
%     num2str(numel_con_LSF_1_LE_0(end)), ' made it to completion'])
% axis([0, Tstore(end-100+1), 0, 40])
% 
% %Number of elements - con LSF_1_LE_p2
% subplot(7,1,2)
% plot(Tstore, numel_con_LSF_1_LE_p2,'k-','LineWidth',2)
% xlabel(['Time (s), con - LSF: 2, (n = ',...
%     num2str(NumOfconLSF_1_LE_p2_Files),') at start, ',...
%     num2str(numel_con_LSF_1_LE_p2(end)), ' made it to completion'])
% axis([0, Tstore(end-100+1), 0, 40])
% 
% %Number of elements - con LSF_1_LE_p4
% subplot(7,1,3)
% plot(Tstore, numel_con_LSF_1_LE_p4,'k-','LineWidth',2)
% xlabel(['Time (s), con - LSF: 3, (n = ',...
%     num2str(NumOfconLSF_1_LE_p4_Files),') at start, ',...
%     num2str(numel_con_LSF_1_LE_p4(end)), ' made it to completion'])
% axis([0, Tstore(end-100+1), 0, 40])
% 
% %Number of elements - con LSF_1_LE_p6
% subplot(7,1,4)
% plot(Tstore, numel_con_LSF_1_LE_p6,'k-','LineWidth',2)
% xlabel(['Time (s), con - LSF: 4, (n = ',...
%     num2str(NumOfconLSF_1_LE_p6_Files),') at start, ',...
%     num2str(numel_con_LSF_1_LE_p6(end)), ' made it to completion'])
% axis([0, Tstore(end-100+1), 0, 40])
% 
% %Number of elements - con LSF_1_LE_p8
% subplot(7,1,5)
% plot(Tstore, numel_con_LSF_1_LE_p8,'k-','LineWidth',2)
% xlabel(['Time (s), con - LSF: 5, (n = ',...
%     num2str(NumOfconLSF_1_LE_p8_Files),') at start, ',...
%     num2str(numel_con_LSF_1_LE_p8(end)), ' made it to completion'])
% axis([0, Tstore(end-100+1), 0, 40])
% 
% %Number of elements - con LSF_1_LE_1
% subplot(7,1,6)
% plot(Tstore, numel_con_LSF_1_LE_1,'k-','LineWidth',2)
% xlabel(['Time (s), con - LSF: 6, (n = ',...
%     num2str(NumOfconLSF_1_LE_1_Files),') at start, ',...
%     num2str(numel_con_LSF_1_LE_1(end)), ' made it to completion'])
% axis([0, Tstore(end-100+1), 0, 40])
% 
% %Number of elements - con LSF_1_LE_2
% subplot(7,1,7)
% plot(Tstore, numel_con_LSF_1_LE_2,'k-','LineWidth',2)
% xlabel(['Time (s), con - LSF: 7, (n = ',...
%     num2str(NumOfconLSF_1_LE_2_Files),') at start, ',...
%     num2str(numel_con_LSF_1_LE_2(end)), ' made it to completion'])
% axis([0, Tstore(end-100+1), 0, 40])
% 
% saveas(gcf,'../Figures_MPC/LSF_1/Fig5_MPC_numel_LSF_1_v12h.fig');
% saveas(gcf,'../Figures_MPC/LSF_1/Fig5_MPC_numel_LSF_1_v12h.jpg');
% 
% disp('Figure 5 done')

%% Cost function violin plots

fig6 = figure(6);
%Cost - con LSF_1_LE_0
subplot(1,7,1)
distributionPlot(mean_cost_con_LSF_1_LE_0,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_0_Files),')'])
ylabel('Mean of cost function value')
% axis([0, 2, cost_min, cost_max])
% set(gca,'ytick',0:2e6:13e6);

%Cost - con LSF_1_LE_p2
subplot(1,7,2)
distributionPlot(mean_cost_con_LSF_1_LE_p2,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.2'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p2_Files),')'])
% axis([0, 2, cost_min, cost_max])
% set(gca,'ytick',0:2e6:13e6);

%Cost - con LSF_1_LE_p4
subplot(1,7,3)
distributionPlot(mean_cost_con_LSF_1_LE_p4,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.4'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p4_Files),')'])
% axis([0, 2, cost_min, cost_max])
% set(gca,'ytick',0:2e6:13e6);

%Cost - con LSF_1_LE_p6
subplot(1,7,4)
distributionPlot(mean_cost_con_LSF_1_LE_p6,'addSpread',true)
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p6_Files),')'])
title({['LSF: ',num2str(LSF_1)];'LE: 0.6'})
% axis([0, 2, cost_min, cost_max])
% set(gca,'ytick',0:2e11:2.5e11);

%Cost - con LSF_1_LE_p8
subplot(1,7,5)
distributionPlot(mean_cost_con_LSF_1_LE_p8,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 0.8'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_p8_Files),')'])
% axis([0, 2, cost_min, cost_max])
% set(gca,'ytick',0:2e11:2.5e11);

%Cost - con LSF_1_LE_1
subplot(1,7,6)
distributionPlot(mean_cost_con_LSF_1_LE_1,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 1'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_1_Files),')'])
% axis([0, 2, cost_min, cost_max])
% set(gca,'ytick',0:2e11:2.5e11);

%Cost - con LSF_1_LE_2
subplot(1,7,7)
distributionPlot(mean_cost_con_LSF_1_LE_2,'addSpread',true)
title({['LSF: ',num2str(LSF_1)];'LE: 2'})
xlabel(['(n=',num2str(NumOfconLSF_1_LE_2_Files),')'])
% axis([0, 2, cost_min, cost_max])
% set(gca,'ytick',0:2e11:2.5e11);

saveas(gcf,'../Figures_MPC/LSF_1/Fig6_MPC_CostViolin_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig6_MPC_CostViolin_LSF_1_v12h.jpg');

set(fig6, 'Position', [120,120,850,525])

disp('Figure 6 done') %p-value = 5.323e-11

%% log plot of cost

figure(7);
semilogy(Tstore(1:timestep:(end-100)),mean_impcost_con_LSF_1_LE_0,'k','Linewidth',2)
hold on;
semilogy(Tstore(1:timestep:(end-100)),mean_impcost_con_LSF_1_LE_p2,'Linewidth',2)
hold on;
semilogy(Tstore(1:timestep:(end-100)),mean_impcost_con_LSF_1_LE_p4,'Linewidth',2)
hold on;
semilogy(Tstore(1:timestep:(end-100)),mean_impcost_con_LSF_1_LE_p6,'Linewidth',2)
hold on;
semilogy(Tstore(1:timestep:(end-100)),mean_impcost_con_LSF_1_LE_p8,'Linewidth',2)
hold on;
semilogy(Tstore(1:timestep:(end-100)),mean_impcost_con_LSF_1_LE_1,'Linewidth',2)
hold on;
semilogy(Tstore(1:timestep:(end-100)),mean_impcost_con_LSF_1_LE_2,'Linewidth',2)
hold on;
xlabel('Time (seconds)'); ylabel('Cost function value')
legend(['LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOfconLSF_1_LE_0_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 0.2 - ', num2str(NumOfconLSF_1_LE_p2_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 0.4 - ', num2str(NumOfconLSF_1_LE_p4_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 0.6 - ', num2str(NumOfconLSF_1_LE_p6_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 0.8 - ', num2str(NumOfconLSF_1_LE_p8_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 1 - ', num2str(NumOfconLSF_1_LE_1_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 2 - ', num2str(NumOfconLSF_1_LE_2_Files)])

saveas(gcf,'../Figures_MPC/LSF_1/Fig7_MPC_LogCost_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig7_MPC_LogCost_LSF_1_v12h.jpg');

disp('Figure 7 done') 

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])