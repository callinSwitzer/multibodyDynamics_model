%8/8/17 Script developed to curate the data into a form more usable by
%error_12h
%8/16/18 Script modified for LengthScaleFactor (LSF)
%10/17/18 Script modified for size extension (LE)
%1/11/19 Script modified for Model Predictive Control (MPC)
    %1/16/19 Script corrected previous indexing error.

%v13b is for horizontal aggressive maneuver
%v13c is for aggressive maneuver
%v12h is for sum of prime number sines

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% List the directory files and assign the number of each file type
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

NumOfconLSF_1_LE_0_Files = numel(listOfconLSF_1_LE_0_Files);
NumOfconLSF_1_LE_p2_Files = numel(listOfconLSF_1_LE_p2_Files);
NumOfconLSF_1_LE_p4_Files = numel(listOfconLSF_1_LE_p4_Files);
NumOfconLSF_1_LE_p6_Files = numel(listOfconLSF_1_LE_p6_Files);
NumOfconLSF_1_LE_p8_Files = numel(listOfconLSF_1_LE_p8_Files);
NumOfconLSF_1_LE_1_Files = numel(listOfconLSF_1_LE_1_Files);
NumOfconLSF_1_LE_2_Files = numel(listOfconLSF_1_LE_2_Files);

%% Assign the imported numbers to internal variables
disp('Assign internal variables')
hws = 500; %Half wingstrokes
PartPath = 2000; %Partial paths (i.e. four partial paths per half wing stroke)
timestep = numel(Tstore(1:(end-100)))/hws;
tderiv = Tstore(2);

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

%% Column bounds for time signal (because of index offset)

bound_a = 2:25:numel(Tstore);
bound_b = 26:25:numel(Tstore);

%% Import LSF_1, LE_0
disp(['Import ',num2str(NumOfconLSF_1_LE_0_Files),' con LSF_1, LE_0 files'])
%State variables of all the successful simulations (i.e. least cost)
Winstore_x_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_y_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_theta_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_phi_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_xdot_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_ydot_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_thetadot_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_phidot_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));
Winstore_beta_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, (hws*timestep+1));

ICsave_am_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files,12*hws);

%Prescribe the first value before the loop.
%Note: x, y, theta dot, and phi dot are all zero initially.
Winstore_theta_con_LSF_1_LE_0(:,1) = pi/4;
Winstore_phi_con_LSF_1_LE_0(:,1) = pi/4 + pi;
Winstore_xdot_con_LSF_1_LE_0(:,1) = 1e-4;
Winstore_ydot_con_LSF_1_LE_0(:,1) = 1e-4;

%State variables
for i = 1:NumOfconLSF_1_LE_0_Files
        
    load(['../SimData_MPC/LSF_1/Winstore_',num2str(i),'_hws_sp_con_LSF_1_LE_0.mat']);
    %State variable extraction and placement
    for j = 1:PartPath
        
        if isfield(Winstore{j,1},['PartPath',num2str(j)]) == 1
            break
        end
        
        ColStartIC = 1+((j-1)*12);
        ColEndIC = ((j-1)*12)+8;
        ColStartVals = 9+((j-1)*12);
        ColEndVals = ((j-1)*12)+11;
        ColCost = 12*j;
        
Winstore_x_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),1)';
        Winstore_y_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),2)';
        Winstore_theta_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),3)'; %in RADIANS
        Winstore_phi_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),4)'; %in RADIANS
        Winstore_xdot_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),5)';
        Winstore_ydot_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),6)';
        Winstore_thetadot_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),7)'; %in RADIANS/second
        Winstore_phidot_con_LSF_1_LE_0(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),8)'; %in RADIANS/second
        
        ICsave_am_con_LSF_1_LE_0(i,ColStartIC:ColEndIC) = Winstore{j,1}.NewICs'; %NewICs
        ICsave_am_con_LSF_1_LE_0(i,ColStartVals:ColEndVals) = Winstore{j,1}.check; %F, alpha, tau
        ICsave_am_con_LSF_1_LE_0(i,ColCost) = Winstore{j,1}.cost; %Cost
    end
end

%% Import LSF_1, LE_p2
disp(['Import ',num2str(NumOfconLSF_1_LE_p2_Files),' con LSF_1, LE_p2 files'])
%State variables of all the successful simulations (i.e. least cost)
Winstore_x_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_y_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_theta_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_phi_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_xdot_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_ydot_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_thetadot_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_phidot_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));
Winstore_beta_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, (hws*timestep+1));

ICsave_am_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files,12*hws);

%Prescribe the first value before the loop.
%Note: x, y, theta dot, and phi dot are all zero initially.
Winstore_theta_con_LSF_1_LE_p2(:,1) = pi/4;
Winstore_phi_con_LSF_1_LE_p2(:,1) = pi/4 + pi;
Winstore_xdot_con_LSF_1_LE_p2(:,1) = 1e-4;
Winstore_ydot_con_LSF_1_LE_p2(:,1) = 1e-4;

%State variables
for i = 1:NumOfconLSF_1_LE_p2_Files
        
    load(['../SimData_MPC/LSF_1/Winstore_',num2str(i),'_hws_sp_con_LSF_1_LE_p2.mat']);
    %State variable extraction and placement
    for j = 1:PartPath
        
        if isfield(Winstore{j,1},['PartPath',num2str(j)]) == 1
            break
        end
        
        ColStartIC = 1+((j-1)*12);
        ColEndIC = ((j-1)*12)+8;
        ColStartVals = 9+((j-1)*12);
        ColEndVals = ((j-1)*12)+11;
        ColCost = 12*j;
        
        Winstore_x_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),1)';
        Winstore_y_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),2)';
        Winstore_theta_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),3)'; %in RADIANS
        Winstore_phi_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),4)'; %in RADIANS
        Winstore_xdot_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),5)';
        Winstore_ydot_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),6)';
        Winstore_thetadot_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),7)'; %in RADIANS/second
        Winstore_phidot_con_LSF_1_LE_p2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),8)'; %in RADIANS/second    
        
        ICsave_am_con_LSF_1_LE_p2(i,ColStartIC:ColEndIC) = Winstore{j,1}.NewICs'; %NewICs
        ICsave_am_con_LSF_1_LE_p2(i,ColStartVals:ColEndVals) = Winstore{j,1}.check; %F, alpha, tau
        ICsave_am_con_LSF_1_LE_p2(i,ColCost) = Winstore{j,1}.cost; %Cost
    end
end

%% Import LSF_1, LE_p4
disp(['Import ',num2str(NumOfconLSF_1_LE_p4_Files),' con LSF_1, LE_p4 files'])
%State variables of all the successful simulations (i.e. least cost)
Winstore_x_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_y_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_theta_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_phi_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_xdot_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_ydot_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_thetadot_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_phidot_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));
Winstore_beta_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, (hws*timestep+1));

ICsave_am_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files,12*hws);

%Prescribe the first value before the loop.
%Note: x, y, theta dot, and phi dot are all zero initially.
Winstore_theta_con_LSF_1_LE_p4(:,1) = pi/4;
Winstore_phi_con_LSF_1_LE_p4(:,1) = pi/4 + pi;
Winstore_xdot_con_LSF_1_LE_p4(:,1) = 1e-4;
Winstore_ydot_con_LSF_1_LE_p4(:,1) = 1e-4;

%State variables
for i = 1:NumOfconLSF_1_LE_p4_Files
        
    load(['../SimData_MPC/LSF_1/Winstore_',num2str(i),'_hws_sp_con_LSF_1_LE_p4.mat']);
    %State variable extraction and placement
    for j = 1:PartPath
        
        if isfield(Winstore{j,1},['PartPath',num2str(j)]) == 1
            break
        end
        
        ColStartIC = 1+((j-1)*12);
        ColEndIC = ((j-1)*12)+8;
        ColStartVals = 9+((j-1)*12);
        ColEndVals = ((j-1)*12)+11;
        ColCost = 12*j;
        
        Winstore_x_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),1)';
        Winstore_y_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),2)';
        Winstore_theta_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),3)'; %in RADIANS
        Winstore_phi_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),4)'; %in RADIANS
        Winstore_xdot_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),5)';
        Winstore_ydot_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),6)';
        Winstore_thetadot_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),7)'; %in RADIANS/second
        Winstore_phidot_con_LSF_1_LE_p4(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),8)'; %in RADIANS/second    
        
        ICsave_am_con_LSF_1_LE_p4(i,ColStartIC:ColEndIC) = Winstore{j,1}.NewICs'; %NewICs
        ICsave_am_con_LSF_1_LE_p4(i,ColStartVals:ColEndVals) = Winstore{j,1}.check; %F, alpha, tau
        ICsave_am_con_LSF_1_LE_p4(i,ColCost) = Winstore{j,1}.cost; %Cost
    end
end

%% Import LSF_1, LE_p6
disp(['Import ',num2str(NumOfconLSF_1_LE_p6_Files),' con LSF_1, LE_p6 files'])
%State variables of all the successful simulations (i.e. least cost)
Winstore_x_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_y_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_theta_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_phi_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_xdot_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_ydot_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_thetadot_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_phidot_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));
Winstore_beta_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, (hws*timestep+1));

ICsave_am_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files,12*hws);

%Prescribe the first value before the loop.
%Note: x, y, theta dot, and phi dot are all zero initially.
Winstore_theta_con_LSF_1_LE_p6(:,1) = pi/4;
Winstore_phi_con_LSF_1_LE_p6(:,1) = pi/4 + pi;
Winstore_xdot_con_LSF_1_LE_p6(:,1) = 1e-4;
Winstore_ydot_con_LSF_1_LE_p6(:,1) = 1e-4;

%State variables
for i = 1:NumOfconLSF_1_LE_p6_Files
        
    load(['../SimData_MPC/LSF_1/Winstore_',num2str(i),'_hws_sp_con_LSF_1_LE_p6.mat']);
    %State variable extraction and placement
    for j = 1:PartPath
        
        if isfield(Winstore{j,1},['PartPath',num2str(j)]) == 1
            break
        end
        
        ColStartIC = 1+((j-1)*12);
        ColEndIC = ((j-1)*12)+8;
        ColStartVals = 9+((j-1)*12);
        ColEndVals = ((j-1)*12)+11;
        ColCost = 12*j;
        
        Winstore_x_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),1)';
        Winstore_y_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),2)';
        Winstore_theta_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),3)'; %in RADIANS
        Winstore_phi_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),4)'; %in RADIANS
        Winstore_xdot_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),5)';
        Winstore_ydot_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),6)';
        Winstore_thetadot_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),7)'; %in RADIANS/second
        Winstore_phidot_con_LSF_1_LE_p6(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),8)'; %in RADIANS/second    
        
        ICsave_am_con_LSF_1_LE_p6(i,ColStartIC:ColEndIC) = Winstore{j,1}.NewICs'; %NewICs
        ICsave_am_con_LSF_1_LE_p6(i,ColStartVals:ColEndVals) = Winstore{j,1}.check; %F, alpha, tau
        ICsave_am_con_LSF_1_LE_p6(i,ColCost) = Winstore{j,1}.cost; %Cost
    end
end

%% Import LSF_1, LE_p8
disp(['Import ',num2str(NumOfconLSF_1_LE_p8_Files),' con LSF_1, LE_p8 files'])
%State variables of all the successful simulations (i.e. least cost)
Winstore_x_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_y_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_theta_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_phi_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_xdot_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_ydot_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_thetadot_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_phidot_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));
Winstore_beta_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, (hws*timestep+1));

ICsave_am_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files,12*hws);

%Prescribe the first value before the loop.
%Note: x, y, theta dot, and phi dot are all zero initially.
Winstore_theta_con_LSF_1_LE_p8(:,1) = pi/4;
Winstore_phi_con_LSF_1_LE_p8(:,1) = pi/4 + pi;
Winstore_xdot_con_LSF_1_LE_p8(:,1) = 1e-4;
Winstore_ydot_con_LSF_1_LE_p8(:,1) = 1e-4;

%State variables
for i = 1:NumOfconLSF_1_LE_p8_Files
        
    load(['../SimData_MPC/LSF_1/Winstore_',num2str(i),'_hws_sp_con_LSF_1_LE_p8.mat']);
    %State variable extraction and placement
    for j = 1:PartPath
        
        if isfield(Winstore{j,1},['PartPath',num2str(j)]) == 1
            break
        end
        
        ColStartIC = 1+((j-1)*12);
        ColEndIC = ((j-1)*12)+8;
        ColStartVals = 9+((j-1)*12);
        ColEndVals = ((j-1)*12)+11;
        ColCost = 12*j;
        
        Winstore_x_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),1)';
        Winstore_y_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),2)';
        Winstore_theta_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),3)'; %in RADIANS
        Winstore_phi_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),4)'; %in RADIANS
        Winstore_xdot_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),5)';
        Winstore_ydot_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),6)';
        Winstore_thetadot_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),7)'; %in RADIANS/second
        Winstore_phidot_con_LSF_1_LE_p8(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),8)'; %in RADIANS/second    
        
        ICsave_am_con_LSF_1_LE_p8(i,ColStartIC:ColEndIC) = Winstore{j,1}.NewICs'; %NewICs
        ICsave_am_con_LSF_1_LE_p8(i,ColStartVals:ColEndVals) = Winstore{j,1}.check; %F, alpha, tau
        ICsave_am_con_LSF_1_LE_p8(i,ColCost) = Winstore{j,1}.cost; %Cost
    end
end

%% Import LSF_1, LE_1
disp(['Import ',num2str(NumOfconLSF_1_LE_1_Files),' con LSF_1, LE_1 files'])
%State variables of all the successful simulations (i.e. least cost)
Winstore_x_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_y_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_theta_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_phi_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_xdot_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_ydot_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_thetadot_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_phidot_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));
Winstore_beta_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, (hws*timestep+1));

ICsave_am_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files,12*hws);

%Prescribe the first value before the loop.
%Note: x, y, theta dot, and phi dot are all zero initially.
Winstore_theta_con_LSF_1_LE_1(:,1) = pi/4;
Winstore_phi_con_LSF_1_LE_1(:,1) = pi/4 + pi;
Winstore_xdot_con_LSF_1_LE_1(:,1) = 1e-4;
Winstore_ydot_con_LSF_1_LE_1(:,1) = 1e-4;

%State variables
for i = 1:NumOfconLSF_1_LE_1_Files
        
    load(['../SimData_MPC/LSF_1/Winstore_',num2str(i),'_hws_sp_con_LSF_1_LE_1.mat']);
    %State variable extraction and placement
    for j = 1:PartPath
        
        if isfield(Winstore{j,1},['PartPath',num2str(j)]) == 1
            break
        end
        
        ColStart = 1+((j-1)*(100*hws/PartPath));
        ColEnd = ((j-1)*(100*hws/PartPath))+(100*hws/PartPath);
        
        ColStartIC = 1+((j-1)*12);
        ColEndIC = ((j-1)*12)+8;
        ColStartVals = 9+((j-1)*12);
        ColEndVals = ((j-1)*12)+11;
        ColCost = 12*j;
        
        Winstore_x_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),1)';
        Winstore_y_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),2)';
        Winstore_theta_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),3)'; %in RADIANS
        Winstore_phi_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),4)'; %in RADIANS
        Winstore_xdot_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),5)';
        Winstore_ydot_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),6)';
        Winstore_thetadot_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),7)'; %in RADIANS/second
        Winstore_phidot_con_LSF_1_LE_1(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),8)'; %in RADIANS/second   
        
        ICsave_am_con_LSF_1_LE_1(i,ColStartIC:ColEndIC) = Winstore{j,1}.NewICs'; %NewICs
        ICsave_am_con_LSF_1_LE_1(i,ColStartVals:ColEndVals) = Winstore{j,1}.check; %F, alpha, tau
        ICsave_am_con_LSF_1_LE_1(i,ColCost) = Winstore{j,1}.cost; %Cost
    end
end

%% Import LSF_1, LE_2
disp(['Import ',num2str(NumOfconLSF_1_LE_2_Files),' con LSF_1, LE_2 files'])
%State variables of all the successful simulations (i.e. least cost)
Winstore_x_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_y_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_theta_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_phi_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_xdot_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_ydot_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_thetadot_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_phidot_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));
Winstore_beta_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, (hws*timestep+1));

ICsave_am_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files,12*hws);

%Prescribe the first value before the loop.
%Note: x, y, theta dot, and phi dot are all zero initially.
Winstore_theta_con_LSF_1_LE_2(:,1) = pi/4;
Winstore_phi_con_LSF_1_LE_2(:,1) = pi/4 + pi;
Winstore_xdot_con_LSF_1_LE_2(:,1) = 1e-4;
Winstore_ydot_con_LSF_1_LE_2(:,1) = 1e-4;

%State variables
for i = 1:NumOfconLSF_1_LE_2_Files
        
    load(['../SimData_MPC/LSF_1/Winstore_',num2str(i),'_hws_sp_con_LSF_1_LE_2.mat']);
    %State variable extraction and placement
    for j = 1:PartPath
        
        if isfield(Winstore{j,1},['PartPath',num2str(j)]) == 1
            break
        end
        
        ColStartIC = 1+((j-1)*12);
        ColEndIC = ((j-1)*12)+8;
        ColStartVals = 9+((j-1)*12);
        ColEndVals = ((j-1)*12)+11;
        ColCost = 12*j;
        
        Winstore_x_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),1)';
        Winstore_y_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),2)';
        Winstore_theta_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),3)'; %in RADIANS
        Winstore_phi_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),4)'; %in RADIANS
        Winstore_xdot_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),5)';
        Winstore_ydot_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),6)';
        Winstore_thetadot_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),7)'; %in RADIANS/second
        Winstore_phidot_con_LSF_1_LE_2(i,bound_a(j):bound_b(j)) =...
            Winstore{j,1}.bigQ(2:(1+100*hws/PartPath),8)'; %in RADIANS/second  
        
        ICsave_am_con_LSF_1_LE_2(i,ColStartIC:ColEndIC) = Winstore{j,1}.NewICs'; %NewICs
        ICsave_am_con_LSF_1_LE_2(i,ColStartVals:ColEndVals) = Winstore{j,1}.check; %F, alpha, tau
        ICsave_am_con_LSF_1_LE_2(i,ColCost) = Winstore{j,1}.cost; %Cost
    end
end

%% To ignore the zeroes that skew our values
disp('Removing unnecessary zeros')

Winstore_x_con_LSF_1_LE_0(Winstore_x_con_LSF_1_LE_0==0) = NaN;
Winstore_y_con_LSF_1_LE_0(Winstore_y_con_LSF_1_LE_0==0) = NaN;
Winstore_theta_con_LSF_1_LE_0(Winstore_theta_con_LSF_1_LE_0==0) = NaN;
Winstore_phi_con_LSF_1_LE_0(Winstore_phi_con_LSF_1_LE_0==0) = NaN;
Winstore_xdot_con_LSF_1_LE_0(Winstore_xdot_con_LSF_1_LE_0==0) = NaN;
Winstore_ydot_con_LSF_1_LE_0(Winstore_ydot_con_LSF_1_LE_0==0) = NaN;
Winstore_thetadot_con_LSF_1_LE_0(Winstore_thetadot_con_LSF_1_LE_0==0) = NaN;
Winstore_phidot_con_LSF_1_LE_0(Winstore_phidot_con_LSF_1_LE_0==0) = NaN;

Winstore_x_con_LSF_1_LE_p2(Winstore_x_con_LSF_1_LE_p2==0) = NaN;
Winstore_y_con_LSF_1_LE_p2(Winstore_y_con_LSF_1_LE_p2==0) = NaN;
Winstore_theta_con_LSF_1_LE_p2(Winstore_theta_con_LSF_1_LE_p2==0) = NaN;
Winstore_phi_con_LSF_1_LE_p2(Winstore_phi_con_LSF_1_LE_p2==0) = NaN;
Winstore_xdot_con_LSF_1_LE_p2(Winstore_xdot_con_LSF_1_LE_p2==0) = NaN;
Winstore_ydot_con_LSF_1_LE_p2(Winstore_ydot_con_LSF_1_LE_p2==0) = NaN;
Winstore_thetadot_con_LSF_1_LE_p2(Winstore_thetadot_con_LSF_1_LE_p2==0) = NaN;
Winstore_phidot_con_LSF_1_LE_p2(Winstore_phidot_con_LSF_1_LE_p2==0) = NaN;

Winstore_x_con_LSF_1_LE_p4(Winstore_x_con_LSF_1_LE_p4==0) = NaN;
Winstore_y_con_LSF_1_LE_p4(Winstore_y_con_LSF_1_LE_p4==0) = NaN;
Winstore_theta_con_LSF_1_LE_p4(Winstore_theta_con_LSF_1_LE_p4==0) = NaN;
Winstore_phi_con_LSF_1_LE_p4(Winstore_phi_con_LSF_1_LE_p4==0) = NaN;
Winstore_xdot_con_LSF_1_LE_p4(Winstore_xdot_con_LSF_1_LE_p4==0) = NaN;
Winstore_ydot_con_LSF_1_LE_p4(Winstore_ydot_con_LSF_1_LE_p4==0) = NaN;
Winstore_thetadot_con_LSF_1_LE_p4(Winstore_thetadot_con_LSF_1_LE_p4==0) = NaN;
Winstore_phidot_con_LSF_1_LE_p4(Winstore_phidot_con_LSF_1_LE_p4==0) = NaN;

Winstore_x_con_LSF_1_LE_p6(Winstore_x_con_LSF_1_LE_p6==0) = NaN;
Winstore_y_con_LSF_1_LE_p6(Winstore_y_con_LSF_1_LE_p6==0) = NaN;
Winstore_theta_con_LSF_1_LE_p6(Winstore_theta_con_LSF_1_LE_p6==0) = NaN;
Winstore_phi_con_LSF_1_LE_p6(Winstore_phi_con_LSF_1_LE_p6==0) = NaN;
Winstore_xdot_con_LSF_1_LE_p6(Winstore_xdot_con_LSF_1_LE_p6==0) = NaN;
Winstore_ydot_con_LSF_1_LE_p6(Winstore_ydot_con_LSF_1_LE_p6==0) = NaN;
Winstore_thetadot_con_LSF_1_LE_p6(Winstore_thetadot_con_LSF_1_LE_p6==0) = NaN;
Winstore_phidot_con_LSF_1_LE_p6(Winstore_phidot_con_LSF_1_LE_p6==0) = NaN;

Winstore_x_con_LSF_1_LE_p8(Winstore_x_con_LSF_1_LE_p8==0) = NaN;
Winstore_y_con_LSF_1_LE_p8(Winstore_y_con_LSF_1_LE_p8==0) = NaN;
Winstore_theta_con_LSF_1_LE_p8(Winstore_theta_con_LSF_1_LE_p8==0) = NaN;
Winstore_phi_con_LSF_1_LE_p8(Winstore_phi_con_LSF_1_LE_p8==0) = NaN;
Winstore_xdot_con_LSF_1_LE_p8(Winstore_xdot_con_LSF_1_LE_p8==0) = NaN;
Winstore_ydot_con_LSF_1_LE_p8(Winstore_ydot_con_LSF_1_LE_p8==0) = NaN;
Winstore_thetadot_con_LSF_1_LE_p8(Winstore_thetadot_con_LSF_1_LE_p8==0) = NaN;
Winstore_phidot_con_LSF_1_LE_p8(Winstore_phidot_con_LSF_1_LE_p8==0) = NaN;

Winstore_x_con_LSF_1_LE_1(Winstore_x_con_LSF_1_LE_1==0) = NaN;
Winstore_y_con_LSF_1_LE_1(Winstore_y_con_LSF_1_LE_1==0) = NaN;
Winstore_theta_con_LSF_1_LE_1(Winstore_theta_con_LSF_1_LE_1==0) = NaN;
Winstore_phi_con_LSF_1_LE_1(Winstore_phi_con_LSF_1_LE_1==0) = NaN;
Winstore_xdot_con_LSF_1_LE_1(Winstore_xdot_con_LSF_1_LE_1==0) = NaN;
Winstore_ydot_con_LSF_1_LE_1(Winstore_ydot_con_LSF_1_LE_1==0) = NaN;
Winstore_thetadot_con_LSF_1_LE_1(Winstore_thetadot_con_LSF_1_LE_1==0) = NaN;
Winstore_phidot_con_LSF_1_LE_1(Winstore_phidot_con_LSF_1_LE_1==0) = NaN;

Winstore_x_con_LSF_1_LE_2(Winstore_x_con_LSF_1_LE_2==0) = NaN;
Winstore_y_con_LSF_1_LE_2(Winstore_y_con_LSF_1_LE_2==0) = NaN;
Winstore_theta_con_LSF_1_LE_2(Winstore_theta_con_LSF_1_LE_2==0) = NaN;
Winstore_phi_con_LSF_1_LE_2(Winstore_phi_con_LSF_1_LE_2==0) = NaN;
Winstore_xdot_con_LSF_1_LE_2(Winstore_xdot_con_LSF_1_LE_2==0) = NaN;
Winstore_ydot_con_LSF_1_LE_2(Winstore_ydot_con_LSF_1_LE_2==0) = NaN;
Winstore_thetadot_con_LSF_1_LE_2(Winstore_thetadot_con_LSF_1_LE_2==0) = NaN;
Winstore_phidot_con_LSF_1_LE_2(Winstore_phidot_con_LSF_1_LE_2==0) = NaN;

%Reassign the initial condition zeros for x, y, thetadot and phidot
Winstore_x_con_LSF_1_LE_0(:,1) = 0;
Winstore_y_con_LSF_1_LE_0(:,1) = 0;
Winstore_thetadot_con_LSF_1_LE_0(:,1) = 0;
Winstore_phidot_con_LSF_1_LE_0(:,1) = 0;

Winstore_x_con_LSF_1_LE_p2(:,1) = 0;
Winstore_y_con_LSF_1_LE_p2(:,1) = 0;
Winstore_thetadot_con_LSF_1_LE_p2(:,1) = 0;
Winstore_phidot_con_LSF_1_LE_p2(:,1) = 0;

Winstore_x_con_LSF_1_LE_p4(:,1) = 0;
Winstore_y_con_LSF_1_LE_p4(:,1) = 0;
Winstore_thetadot_con_LSF_1_LE_p4(:,1) = 0;
Winstore_phidot_con_LSF_1_LE_p4(:,1) = 0;

Winstore_x_con_LSF_1_LE_p6(:,1) = 0;
Winstore_y_con_LSF_1_LE_p6(:,1) = 0;
Winstore_thetadot_con_LSF_1_LE_p6(:,1) = 0;
Winstore_phidot_con_LSF_1_LE_p6(:,1) = 0;

Winstore_x_con_LSF_1_LE_p8(:,1) = 0;
Winstore_y_con_LSF_1_LE_p8(:,1) = 0;
Winstore_thetadot_con_LSF_1_LE_p8(:,1) = 0;
Winstore_phidot_con_LSF_1_LE_p8(:,1) = 0;

Winstore_x_con_LSF_1_LE_1(:,1) = 0;
Winstore_y_con_LSF_1_LE_1(:,1) = 0;
Winstore_thetadot_con_LSF_1_LE_1(:,1) = 0;
Winstore_phidot_con_LSF_1_LE_1(:,1) = 0;

Winstore_x_con_LSF_1_LE_2(:,1) = 0;
Winstore_y_con_LSF_1_LE_2(:,1) = 0;
Winstore_thetadot_con_LSF_1_LE_2(:,1) = 0;
Winstore_phidot_con_LSF_1_LE_2(:,1) = 0;

%% Calculate beta (abdominal flexion) values

Winstore_beta_con_LSF_1_LE_0 = Winstore_phi_con_LSF_1_LE_0 - Winstore_theta_con_LSF_1_LE_0 - pi;
Winstore_beta_con_LSF_1_LE_p2 = Winstore_phi_con_LSF_1_LE_p2 - Winstore_theta_con_LSF_1_LE_p2 - pi;
Winstore_beta_con_LSF_1_LE_p4 = Winstore_phi_con_LSF_1_LE_p4 - Winstore_theta_con_LSF_1_LE_p4 - pi;
Winstore_beta_con_LSF_1_LE_p6 = Winstore_phi_con_LSF_1_LE_p6 - Winstore_theta_con_LSF_1_LE_p6 - pi;
Winstore_beta_con_LSF_1_LE_p8 = Winstore_phi_con_LSF_1_LE_p8 - Winstore_theta_con_LSF_1_LE_p8 - pi;
Winstore_beta_con_LSF_1_LE_1 = Winstore_phi_con_LSF_1_LE_1 - Winstore_theta_con_LSF_1_LE_1 - pi;
Winstore_beta_con_LSF_1_LE_2 = Winstore_phi_con_LSF_1_LE_2 - Winstore_theta_con_LSF_1_LE_2 - pi;

%% Save the Winstore_y data

save('../CuratedData_MPC/LSF_1/Winstore_y_con_LSF_1_LE_0.mat','Winstore_y_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/Winstore_y_con_LSF_1_LE_p2.mat','Winstore_y_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/Winstore_y_con_LSF_1_LE_p4.mat','Winstore_y_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/Winstore_y_con_LSF_1_LE_p6.mat','Winstore_y_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/Winstore_y_con_LSF_1_LE_p8.mat','Winstore_y_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/Winstore_y_con_LSF_1_LE_1.mat','Winstore_y_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/Winstore_y_con_LSF_1_LE_2.mat','Winstore_y_con_LSF_1_LE_2');

disp('Winstore_y data saved')

%% Calculate tracking error
disp('Calculate tracking error')

x_trackerr_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, hws*timestep+1);
y_trackerr_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, hws*timestep+1);
trackingerror_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, hws*timestep+1);
x_trackerr_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, hws*timestep+1);
y_trackerr_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, hws*timestep+1);
trackingerror_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, hws*timestep+1);
x_trackerr_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, hws*timestep+1);
y_trackerr_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, hws*timestep+1);
trackingerror_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, hws*timestep+1);
x_trackerr_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, hws*timestep+1);
y_trackerr_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, hws*timestep+1);
trackingerror_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, hws*timestep+1);
x_trackerr_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, hws*timestep+1);
y_trackerr_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, hws*timestep+1);
trackingerror_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, hws*timestep+1);
x_trackerr_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, hws*timestep+1);
y_trackerr_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, hws*timestep+1);
trackingerror_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, hws*timestep+1);
x_trackerr_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, hws*timestep+1);
y_trackerr_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, hws*timestep+1);
trackingerror_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, hws*timestep+1);

for i = 1:hws*timestep+1
    x_trackerr_con_LSF_1_LE_0(:,i) = Winstore_x_con_LSF_1_LE_0(:,i) - x_g(1,i);
    y_trackerr_con_LSF_1_LE_0(:,i) = Winstore_y_con_LSF_1_LE_0(:,i) - y_g(1,i);
    x_trackerr_con_LSF_1_LE_p2(:,i) = Winstore_x_con_LSF_1_LE_p2(:,i) - x_g(1,i);
    y_trackerr_con_LSF_1_LE_p2(:,i) = Winstore_y_con_LSF_1_LE_p2(:,i) - y_g(1,i);
    x_trackerr_con_LSF_1_LE_p4(:,i) = Winstore_x_con_LSF_1_LE_p4(:,i) - x_g(1,i);
    y_trackerr_con_LSF_1_LE_p4(:,i) = Winstore_y_con_LSF_1_LE_p4(:,i) - y_g(1,i);
    x_trackerr_con_LSF_1_LE_p6(:,i) = Winstore_x_con_LSF_1_LE_p6(:,i) - x_g(1,i);
    y_trackerr_con_LSF_1_LE_p6(:,i) = Winstore_y_con_LSF_1_LE_p6(:,i) - y_g(1,i);
    x_trackerr_con_LSF_1_LE_p8(:,i) = Winstore_x_con_LSF_1_LE_p8(:,i) - x_g(1,i);
    y_trackerr_con_LSF_1_LE_p8(:,i) = Winstore_y_con_LSF_1_LE_p8(:,i) - y_g(1,i);
    x_trackerr_con_LSF_1_LE_1(:,i) = Winstore_x_con_LSF_1_LE_1(:,i) - x_g(1,i);
    y_trackerr_con_LSF_1_LE_1(:,i) = Winstore_y_con_LSF_1_LE_1(:,i) - y_g(1,i);
    x_trackerr_con_LSF_1_LE_2(:,i) = Winstore_x_con_LSF_1_LE_2(:,i) - x_g(1,i);
    y_trackerr_con_LSF_1_LE_2(:,i) = Winstore_y_con_LSF_1_LE_2(:,i) - y_g(1,i);
end

trackingerror_con_LSF_1_LE_0 = sqrt((x_trackerr_con_LSF_1_LE_0).^2 + (y_trackerr_con_LSF_1_LE_0).^2);
trackingerror_con_LSF_1_LE_p2 = sqrt((x_trackerr_con_LSF_1_LE_p2).^2 + (y_trackerr_con_LSF_1_LE_p2).^2);
trackingerror_con_LSF_1_LE_p4 = sqrt((x_trackerr_con_LSF_1_LE_p4).^2 + (y_trackerr_con_LSF_1_LE_p4).^2);
trackingerror_con_LSF_1_LE_p6 = sqrt((x_trackerr_con_LSF_1_LE_p6).^2 + (y_trackerr_con_LSF_1_LE_p6).^2);
trackingerror_con_LSF_1_LE_p8 = sqrt((x_trackerr_con_LSF_1_LE_p8).^2 + (y_trackerr_con_LSF_1_LE_p8).^2);
trackingerror_con_LSF_1_LE_1 = sqrt((x_trackerr_con_LSF_1_LE_1).^2 + (y_trackerr_con_LSF_1_LE_1).^2);
trackingerror_con_LSF_1_LE_2 = sqrt((x_trackerr_con_LSF_1_LE_2).^2 + (y_trackerr_con_LSF_1_LE_2).^2);

%Save the tracking error data
save('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_0.mat','trackingerror_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p2.mat','trackingerror_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p4.mat','trackingerror_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p6.mat','trackingerror_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_p8.mat','trackingerror_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_1.mat','trackingerror_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/trackingerror_con_LSF_1_LE_2.mat','trackingerror_con_LSF_1_LE_2');

%% Calculate mean of tracking error per full run

%tepfr is tracking error per full run
mean_tepfr_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files,1);
mean_tepfr_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files,1);
mean_tepfr_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files,1);
mean_tepfr_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files,1);
mean_tepfr_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files,1);
mean_tepfr_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files,1);
mean_tepfr_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files,1);

for i = 1:NumOfconLSF_1_LE_0_Files
    mean_tepfr_con_LSF_1_LE_0(i,1) = nanmean(trackingerror_con_LSF_1_LE_0(i,:));
end

for i = 1:NumOfconLSF_1_LE_p2_Files
    mean_tepfr_con_LSF_1_LE_p2(i,1) = nanmean(trackingerror_con_LSF_1_LE_p2(i,:));
end

for i = 1:NumOfconLSF_1_LE_p4_Files
    mean_tepfr_con_LSF_1_LE_p4(i,1) = nanmean(trackingerror_con_LSF_1_LE_p4(i,:));
end

for i = 1:NumOfconLSF_1_LE_p6_Files
    mean_tepfr_con_LSF_1_LE_p6(i,1) = nanmean(trackingerror_con_LSF_1_LE_p6(i,:));
end

for i = 1:NumOfconLSF_1_LE_p8_Files
    mean_tepfr_con_LSF_1_LE_p8(i,1) = nanmean(trackingerror_con_LSF_1_LE_p8(i,:));
end

for i = 1:NumOfconLSF_1_LE_1_Files
    mean_tepfr_con_LSF_1_LE_1(i,1) = nanmean(trackingerror_con_LSF_1_LE_1(i,:));
end

for i = 1:NumOfconLSF_1_LE_2_Files
    mean_tepfr_con_LSF_1_LE_2(i,1) = nanmean(trackingerror_con_LSF_1_LE_2(i,:));
end

%Curate for Stats csv file
size_export_con_LSF_1_LE_0(1:numel(mean_tepfr_con_LSF_1_LE_0),1) = {'LSF_1'};
treatment_export_con_LSF_1_LE_0(1:numel(mean_tepfr_con_LSF_1_LE_0),1) = {'con'};
size_export_con_LSF_1_LE_p2(1:numel(mean_tepfr_con_LSF_1_LE_p2),1) = {'LSF_1'};
treatment_export_con_LSF_1_LE_p2(1:numel(mean_tepfr_con_LSF_1_LE_p2),1) = {'con'};
size_export_con_LSF_1_LE_p4(1:numel(mean_tepfr_con_LSF_1_LE_p4),1) = {'LSF_1'};
treatment_export_con_LSF_1_LE_p4(1:numel(mean_tepfr_con_LSF_1_LE_p4),1) = {'con'};
size_export_con_LSF_1_LE_p6(1:numel(mean_tepfr_con_LSF_1_LE_p6),1) = {'LSF_1'};
treatment_export_con_LSF_1_LE_p6(1:numel(mean_tepfr_con_LSF_1_LE_p6),1) = {'con'};
size_export_con_LSF_1_LE_p8(1:numel(mean_tepfr_con_LSF_1_LE_p8),1) = {'LSF_1'};
treatment_export_con_LSF_1_LE_p8(1:numel(mean_tepfr_con_LSF_1_LE_p8),1) = {'con'};
size_export_con_LSF_1_LE_1(1:numel(mean_tepfr_con_LSF_1_LE_1),1) = {'LSF_1'};
treatment_export_con_LSF_1_LE_1(1:numel(mean_tepfr_con_LSF_1_LE_1),1) = {'con'};
size_export_con_LSF_1_LE_2(1:numel(mean_tepfr_con_LSF_1_LE_2),1) = {'LSF_1'};
treatment_export_con_LSF_1_LE_2(1:numel(mean_tepfr_con_LSF_1_LE_2),1) = {'con'};

LE_export_con_LSF_1_LE_0(1:numel(mean_tepfr_con_LSF_1_LE_0),1) = {'LE_0'};
LE_export_con_LSF_1_LE_p2(1:numel(mean_tepfr_con_LSF_1_LE_p2),1) = {'LE_p2'};
LE_export_con_LSF_1_LE_p4(1:numel(mean_tepfr_con_LSF_1_LE_p4),1) = {'LE_p4'};
LE_export_con_LSF_1_LE_p6(1:numel(mean_tepfr_con_LSF_1_LE_p6),1) = {'LE_p6'};
LE_export_con_LSF_1_LE_p8(1:numel(mean_tepfr_con_LSF_1_LE_p8),1) = {'LE_p8'};
LE_export_con_LSF_1_LE_1(1:numel(mean_tepfr_con_LSF_1_LE_1),1) = {'LE_1'};
LE_export_con_LSF_1_LE_2(1:numel(mean_tepfr_con_LSF_1_LE_2),1) = {'LE_2'};


export_size = [size_export_con_LSF_1_LE_0;...
    size_export_con_LSF_1_LE_p2; size_export_con_LSF_1_LE_p4;...
    size_export_con_LSF_1_LE_p6; size_export_con_LSF_1_LE_p8;...
    size_export_con_LSF_1_LE_1; size_export_con_LSF_1_LE_2];

export_treatment = [treatment_export_con_LSF_1_LE_0;...
    treatment_export_con_LSF_1_LE_p2; treatment_export_con_LSF_1_LE_p4;...
    treatment_export_con_LSF_1_LE_p6; treatment_export_con_LSF_1_LE_p8;...
    treatment_export_con_LSF_1_LE_1; treatment_export_con_LSF_1_LE_2];

export_LE = [LE_export_con_LSF_1_LE_0;...
    LE_export_con_LSF_1_LE_p2; LE_export_con_LSF_1_LE_p4;...
    LE_export_con_LSF_1_LE_p6; LE_export_con_LSF_1_LE_p8;...
    LE_export_con_LSF_1_LE_1; LE_export_con_LSF_1_LE_2];

tepfr_export_numbers = [mean_tepfr_con_LSF_1_LE_0;...
    mean_tepfr_con_LSF_1_LE_p2; mean_tepfr_con_LSF_1_LE_p4;...
    mean_tepfr_con_LSF_1_LE_p6; mean_tepfr_con_LSF_1_LE_p8;...
    mean_tepfr_con_LSF_1_LE_1; mean_tepfr_con_LSF_1_LE_2];

tepfr_write = fopen('../StatsStuff_MPC/hws_sp_tepfr_LSF_1_LE.csv','w');
%Write the names to a CSV file
for k = 1:numel(export_size)
    fprintf(tepfr_write,'%f,',tepfr_export_numbers(k,1)); 
                    %The %f is floating point precision number,
                    %the \n indicates a new line
    fprintf(tepfr_write,'%s,', export_size{k,1}); 
                    %The %s is for string formatting,
                    %Make sure the comma is where it is!
    fprintf(tepfr_write,'%s,', export_treatment{k,1}); 
                    %The %s is for string formatting,
                    %Make sure the comma is where it is!
    fprintf(tepfr_write,'%s\n', export_LE{k,1}); 
                    %The %s is for string formatting,
                    %Make sure the comma is where it is!
end
fclose(tepfr_write);

%% Cost function as imported by the earlier values
disp('Cost matrix (as imported)')

cost_import_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files, hws);
cost_import_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files, hws);
cost_import_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files, hws);
cost_import_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files, hws);
cost_import_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files, hws);
cost_import_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files, hws);
cost_import_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files, hws);

for i = 1:hws
    cost_import_con_LSF_1_LE_0(:,i) = ICsave_am_con_LSF_1_LE_0(:,12*i);
    cost_import_con_LSF_1_LE_p2(:,i) = ICsave_am_con_LSF_1_LE_p2(:,12*i);
    cost_import_con_LSF_1_LE_p4(:,i) = ICsave_am_con_LSF_1_LE_p4(:,12*i);
    cost_import_con_LSF_1_LE_p6(:,i) = ICsave_am_con_LSF_1_LE_p6(:,12*i);
    cost_import_con_LSF_1_LE_p8(:,i) = ICsave_am_con_LSF_1_LE_p8(:,12*i);
    cost_import_con_LSF_1_LE_1(:,i) = ICsave_am_con_LSF_1_LE_1(:,12*i);
    cost_import_con_LSF_1_LE_2(:,i) = ICsave_am_con_LSF_1_LE_2(:,12*i);
end

cost_import_con_LSF_1_LE_0(cost_import_con_LSF_1_LE_0==0) = NaN;
cost_import_con_LSF_1_LE_p2(cost_import_con_LSF_1_LE_p2==0) = NaN;
cost_import_con_LSF_1_LE_p4(cost_import_con_LSF_1_LE_p4==0) = NaN;
cost_import_con_LSF_1_LE_p6(cost_import_con_LSF_1_LE_p6==0) = NaN;
cost_import_con_LSF_1_LE_p8(cost_import_con_LSF_1_LE_p8==0) = NaN;
cost_import_con_LSF_1_LE_1(cost_import_con_LSF_1_LE_1==0) = NaN;
cost_import_con_LSF_1_LE_2(cost_import_con_LSF_1_LE_2==0) = NaN;

%% Calculate mean of cost per full run
%Cost
mean_cost_con_LSF_1_LE_0 = zeros(NumOfconLSF_1_LE_0_Files,1);
mean_cost_con_LSF_1_LE_p2 = zeros(NumOfconLSF_1_LE_p2_Files,1);
mean_cost_con_LSF_1_LE_p4 = zeros(NumOfconLSF_1_LE_p4_Files,1);
mean_cost_con_LSF_1_LE_p6 = zeros(NumOfconLSF_1_LE_p6_Files,1);
mean_cost_con_LSF_1_LE_p8 = zeros(NumOfconLSF_1_LE_p8_Files,1);
mean_cost_con_LSF_1_LE_1 = zeros(NumOfconLSF_1_LE_1_Files,1);
mean_cost_con_LSF_1_LE_2 = zeros(NumOfconLSF_1_LE_2_Files,1);

for i = 1:NumOfconLSF_1_LE_0_Files
    mean_cost_con_LSF_1_LE_0(i,1) = nanmean(cost_import_con_LSF_1_LE_0(i,:));
end

for i = 1:NumOfconLSF_1_LE_p2_Files
    mean_cost_con_LSF_1_LE_p2(i,1) = nanmean(cost_import_con_LSF_1_LE_p2(i,:));
end

for i = 1:NumOfconLSF_1_LE_p4_Files
    mean_cost_con_LSF_1_LE_p4(i,1) = nanmean(cost_import_con_LSF_1_LE_p4(i,:));
end

for i = 1:NumOfconLSF_1_LE_p6_Files
    mean_cost_con_LSF_1_LE_p6(i,1) = nanmean(cost_import_con_LSF_1_LE_p6(i,:));
end

for i = 1:NumOfconLSF_1_LE_p8_Files
    mean_cost_con_LSF_1_LE_p8(i,1) = nanmean(cost_import_con_LSF_1_LE_p8(i,:));
end

for i = 1:NumOfconLSF_1_LE_1_Files
    mean_cost_con_LSF_1_LE_1(i,1) = nanmean(cost_import_con_LSF_1_LE_1(i,:));
end

for i = 1:NumOfconLSF_1_LE_2_Files
    mean_cost_con_LSF_1_LE_2(i,1) = nanmean(cost_import_con_LSF_1_LE_2(i,:));
end

cost_export_numbers_LE = [mean_cost_con_LSF_1_LE_0;...
    mean_cost_con_LSF_1_LE_p2; mean_cost_con_LSF_1_LE_p4;...
    mean_cost_con_LSF_1_LE_p6; mean_cost_con_LSF_1_LE_p8;...
    mean_cost_con_LSF_1_LE_1; mean_cost_con_LSF_1_LE_2];

cost_write = fopen('../StatsStuff_MPC/hws_sp_cost_LSF_1_LE.csv','w');
%Write the names to a CSV file
for k = 1:numel(export_size)
    fprintf(cost_write,'%f,',cost_export_numbers_LE(k,1)); 
                    %The %f is floating point precision number,
                    %the \n indicates a new line
    fprintf(cost_write,'%s,', export_size{k,1}); 
                    %The %s is for string formatting,
                    %Make sure the comma is where it is!
    fprintf(cost_write,'%s,', export_treatment{k,1}); 
                    %The %s is for string formatting,
                    %Make sure the comma is where it is!
    fprintf(tepfr_write,'%s\n', export_LE{k,1}); 
                    %The %s is for string formatting,
                    %Make sure the comma is where it is!                
end
fclose(cost_write);

%% Calculate the mean and standard deviation of the variables
disp('Calculate means and st devs')
tic
%LSF_1_LE_0
%Mean of state vars
mean_x_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_y_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_theta_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_phi_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_xdot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_ydot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_thetadot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_phidot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_beta_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
mean_dist_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));

%Standard deviation of state vars
std_x_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_y_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_theta_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_phi_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_xdot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_ydot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_thetadot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_phidot_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_beta_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
std_dist_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));

%LSF_1_LE_p2
%Mean of state vars
mean_x_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_y_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_theta_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_phi_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_xdot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_ydot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_thetadot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_phidot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_beta_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
mean_dist_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));

%Standard deviation of state vars
std_x_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_y_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_theta_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_phi_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_xdot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_ydot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_thetadot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_phidot_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_beta_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
std_dist_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));

%LSF_1_LE_p4
%Mean of state vars
mean_x_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_y_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_theta_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_phi_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_xdot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_ydot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_thetadot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_phidot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_beta_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
mean_dist_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));

%Standard deviation of state vars
std_x_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_y_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_theta_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_phi_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_xdot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_ydot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_thetadot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_phidot_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_beta_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
std_dist_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));

%LSF_1_LE_p6
%Mean of state vars
mean_x_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_y_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_theta_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_phi_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_xdot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_ydot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_thetadot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_phidot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_beta_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
mean_dist_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));

%Standard deviation of state vars
std_x_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_y_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_theta_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_phi_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_xdot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_ydot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_thetadot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_phidot_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_beta_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
std_dist_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));

%LSF_1_LE_p8
%Mean of state vars
mean_x_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_y_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_theta_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_phi_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_xdot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_ydot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_thetadot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_phidot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_beta_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
mean_dist_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));

%Standard deviation of state vars
std_x_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_y_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_theta_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_phi_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_xdot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_ydot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_thetadot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_phidot_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_beta_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
std_dist_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));

%LSF_1_LE_1
%Mean of state vars
mean_x_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_y_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_theta_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_phi_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_xdot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_ydot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_thetadot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_phidot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_beta_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
mean_dist_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));

%Standard deviation of state vars
std_x_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_y_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_theta_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_phi_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_xdot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_ydot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_thetadot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_phidot_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_beta_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
std_dist_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));

%LSF_1_LE_2
%Mean of state vars
mean_x_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_y_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_theta_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_phi_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_xdot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_ydot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_thetadot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_phidot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_beta_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
mean_dist_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));

%Standard deviation of state vars
std_x_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_y_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_theta_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_phi_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_xdot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_ydot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_thetadot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_phidot_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_beta_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));
std_dist_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));

%Calculate number of elements of each full run
numel_con_LSF_1_LE_0 = zeros(1,(hws*timestep+1));
numel_con_LSF_1_LE_p2 = zeros(1,(hws*timestep+1));
numel_con_LSF_1_LE_p4 = zeros(1,(hws*timestep+1));
numel_con_LSF_1_LE_p6 = zeros(1,(hws*timestep+1));
numel_con_LSF_1_LE_p8 = zeros(1,(hws*timestep+1));
numel_con_LSF_1_LE_1 = zeros(1,(hws*timestep+1));
numel_con_LSF_1_LE_2 = zeros(1,(hws*timestep+1));

%Imported cost mean & standard deviations
mean_impcost_con_LSF_1_LE_0 = zeros(1,hws); 
mean_impcost_con_LSF_1_LE_p2 = zeros(1,hws); 
mean_impcost_con_LSF_1_LE_p4 = zeros(1,hws); 
mean_impcost_con_LSF_1_LE_p6 = zeros(1,hws); 
mean_impcost_con_LSF_1_LE_p8 = zeros(1,hws); 
mean_impcost_con_LSF_1_LE_1 = zeros(1,hws); 
mean_impcost_con_LSF_1_LE_2 = zeros(1,hws); 

std_impcost_con_LSF_1_LE_0 = zeros(1,hws);
std_impcost_con_LSF_1_LE_p2 = zeros(1,hws);
std_impcost_con_LSF_1_LE_p4 = zeros(1,hws);
std_impcost_con_LSF_1_LE_p6 = zeros(1,hws);
std_impcost_con_LSF_1_LE_p8 = zeros(1,hws);
std_impcost_con_LSF_1_LE_1 = zeros(1,hws);
std_impcost_con_LSF_1_LE_2 = zeros(1,hws);

for i = 1:(hws*timestep+1)
    mean_x_con_LSF_1_LE_0(1,i) = nanmean(Winstore_x_con_LSF_1_LE_0(:,i));
    mean_y_con_LSF_1_LE_0(1,i) = nanmean(Winstore_y_con_LSF_1_LE_0(:,i));
    mean_theta_con_LSF_1_LE_0(1,i) = nanmean(Winstore_theta_con_LSF_1_LE_0(:,i));
    mean_phi_con_LSF_1_LE_0(1,i) = nanmean(Winstore_phi_con_LSF_1_LE_0(:,i));
    mean_xdot_con_LSF_1_LE_0(1,i) = nanmean(Winstore_xdot_con_LSF_1_LE_0(:,i));
    mean_ydot_con_LSF_1_LE_0(1,i) = nanmean(Winstore_ydot_con_LSF_1_LE_0(:,i));
    mean_thetadot_con_LSF_1_LE_0(1,i) = nanmean(Winstore_thetadot_con_LSF_1_LE_0(:,i));
    mean_phidot_con_LSF_1_LE_0(1,i) = nanmean(Winstore_phidot_con_LSF_1_LE_0(:,i));
    mean_beta_con_LSF_1_LE_0(1,i) = nanmean(Winstore_beta_con_LSF_1_LE_0(:,i));
    mean_dist_con_LSF_1_LE_0(1,i) = nanmean(trackingerror_con_LSF_1_LE_0(:,i));
    std_x_con_LSF_1_LE_0(1,i) = nanstd(Winstore_x_con_LSF_1_LE_0(:,i));
    std_y_con_LSF_1_LE_0(1,i) = nanstd(Winstore_y_con_LSF_1_LE_0(:,i));
    std_theta_con_LSF_1_LE_0(1,i) = nanstd(Winstore_theta_con_LSF_1_LE_0(:,i));
    std_phi_con_LSF_1_LE_0(1,i) = nanstd(Winstore_phi_con_LSF_1_LE_0(:,i));
    std_xdot_con_LSF_1_LE_0(1,i) = nanstd(Winstore_xdot_con_LSF_1_LE_0(:,i));
    std_ydot_con_LSF_1_LE_0(1,i) = nanstd(Winstore_ydot_con_LSF_1_LE_0(:,i));
    std_thetadot_con_LSF_1_LE_0(1,i) = nanstd(Winstore_thetadot_con_LSF_1_LE_0(:,i));
    std_phidot_con_LSF_1_LE_0(1,i) = nanstd(Winstore_phidot_con_LSF_1_LE_0(:,i));
    std_beta_con_LSF_1_LE_0(1,i) = nanstd(Winstore_beta_con_LSF_1_LE_0(:,i));
    std_dist_con_LSF_1_LE_0(1,i) = nanstd(trackingerror_con_LSF_1_LE_0(:,i));
    
    mean_x_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_x_con_LSF_1_LE_p2(:,i));
    mean_y_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_y_con_LSF_1_LE_p2(:,i));
    mean_theta_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_theta_con_LSF_1_LE_p2(:,i));
    mean_phi_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_phi_con_LSF_1_LE_p2(:,i));
    mean_xdot_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_xdot_con_LSF_1_LE_p2(:,i));
    mean_ydot_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_ydot_con_LSF_1_LE_p2(:,i));
    mean_thetadot_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_thetadot_con_LSF_1_LE_p2(:,i));
    mean_phidot_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_phidot_con_LSF_1_LE_p2(:,i));
    mean_beta_con_LSF_1_LE_p2(1,i) = nanmean(Winstore_beta_con_LSF_1_LE_p2(:,i));
    mean_dist_con_LSF_1_LE_p2(1,i) = nanmean(trackingerror_con_LSF_1_LE_p2(:,i));
    std_x_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_x_con_LSF_1_LE_p2(:,i));
    std_y_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_y_con_LSF_1_LE_p2(:,i));
    std_theta_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_theta_con_LSF_1_LE_p2(:,i));
    std_phi_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_phi_con_LSF_1_LE_p2(:,i));
    std_xdot_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_xdot_con_LSF_1_LE_p2(:,i));
    std_ydot_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_ydot_con_LSF_1_LE_p2(:,i));
    std_thetadot_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_thetadot_con_LSF_1_LE_p2(:,i));
    std_phidot_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_phidot_con_LSF_1_LE_p2(:,i));
    std_beta_con_LSF_1_LE_p2(1,i) = nanstd(Winstore_beta_con_LSF_1_LE_p2(:,i));
    std_dist_con_LSF_1_LE_p2(1,i) = nanstd(trackingerror_con_LSF_1_LE_p2(:,i));
    
    mean_x_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_x_con_LSF_1_LE_p4(:,i));
    mean_y_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_y_con_LSF_1_LE_p4(:,i));
    mean_theta_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_theta_con_LSF_1_LE_p4(:,i));
    mean_phi_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_phi_con_LSF_1_LE_p4(:,i));
    mean_xdot_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_xdot_con_LSF_1_LE_p4(:,i));
    mean_ydot_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_ydot_con_LSF_1_LE_p4(:,i));
    mean_thetadot_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_thetadot_con_LSF_1_LE_p4(:,i));
    mean_phidot_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_phidot_con_LSF_1_LE_p4(:,i));
    mean_beta_con_LSF_1_LE_p4(1,i) = nanmean(Winstore_beta_con_LSF_1_LE_p4(:,i));
    mean_dist_con_LSF_1_LE_p4(1,i) = nanmean(trackingerror_con_LSF_1_LE_p4(:,i));
    std_x_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_x_con_LSF_1_LE_p4(:,i));
    std_y_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_y_con_LSF_1_LE_p4(:,i));
    std_theta_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_theta_con_LSF_1_LE_p4(:,i));
    std_phi_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_phi_con_LSF_1_LE_p4(:,i));
    std_xdot_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_xdot_con_LSF_1_LE_p4(:,i));
    std_ydot_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_ydot_con_LSF_1_LE_p4(:,i));
    std_thetadot_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_thetadot_con_LSF_1_LE_p4(:,i));
    std_phidot_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_phidot_con_LSF_1_LE_p4(:,i));
    std_beta_con_LSF_1_LE_p4(1,i) = nanstd(Winstore_beta_con_LSF_1_LE_p4(:,i));
    std_dist_con_LSF_1_LE_p4(1,i) = nanstd(trackingerror_con_LSF_1_LE_p4(:,i));
    
    mean_x_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_x_con_LSF_1_LE_p6(:,i));
    mean_y_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_y_con_LSF_1_LE_p6(:,i));
    mean_theta_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_theta_con_LSF_1_LE_p6(:,i));
    mean_phi_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_phi_con_LSF_1_LE_p6(:,i));
    mean_xdot_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_xdot_con_LSF_1_LE_p6(:,i));
    mean_ydot_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_ydot_con_LSF_1_LE_p6(:,i));
    mean_thetadot_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_thetadot_con_LSF_1_LE_p6(:,i));
    mean_phidot_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_phidot_con_LSF_1_LE_p6(:,i));
    mean_beta_con_LSF_1_LE_p6(1,i) = nanmean(Winstore_beta_con_LSF_1_LE_p6(:,i));
    mean_dist_con_LSF_1_LE_p6(1,i) = nanmean(trackingerror_con_LSF_1_LE_p6(:,i));
    std_x_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_x_con_LSF_1_LE_p6(:,i));
    std_y_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_y_con_LSF_1_LE_p6(:,i));
    std_theta_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_theta_con_LSF_1_LE_p6(:,i));
    std_phi_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_phi_con_LSF_1_LE_p6(:,i));
    std_xdot_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_xdot_con_LSF_1_LE_p6(:,i));
    std_ydot_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_ydot_con_LSF_1_LE_p6(:,i));
    std_thetadot_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_thetadot_con_LSF_1_LE_p6(:,i));
    std_phidot_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_phidot_con_LSF_1_LE_p6(:,i));
    std_beta_con_LSF_1_LE_p6(1,i) = nanstd(Winstore_beta_con_LSF_1_LE_p6(:,i));
    std_dist_con_LSF_1_LE_p6(1,i) = nanstd(trackingerror_con_LSF_1_LE_p6(:,i));
    
    mean_x_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_x_con_LSF_1_LE_p8(:,i));
    mean_y_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_y_con_LSF_1_LE_p8(:,i));
    mean_theta_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_theta_con_LSF_1_LE_p8(:,i));
    mean_phi_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_phi_con_LSF_1_LE_p8(:,i));
    mean_xdot_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_xdot_con_LSF_1_LE_p8(:,i));
    mean_ydot_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_ydot_con_LSF_1_LE_p8(:,i));
    mean_thetadot_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_thetadot_con_LSF_1_LE_p8(:,i));
    mean_phidot_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_phidot_con_LSF_1_LE_p8(:,i));
    mean_beta_con_LSF_1_LE_p8(1,i) = nanmean(Winstore_beta_con_LSF_1_LE_p8(:,i));
    mean_dist_con_LSF_1_LE_p8(1,i) = nanmean(trackingerror_con_LSF_1_LE_p8(:,i));
    std_x_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_x_con_LSF_1_LE_p8(:,i));
    std_y_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_y_con_LSF_1_LE_p8(:,i));
    std_theta_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_theta_con_LSF_1_LE_p8(:,i));
    std_phi_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_phi_con_LSF_1_LE_p8(:,i));
    std_xdot_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_xdot_con_LSF_1_LE_p8(:,i));
    std_ydot_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_ydot_con_LSF_1_LE_p8(:,i));
    std_thetadot_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_thetadot_con_LSF_1_LE_p8(:,i));
    std_phidot_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_phidot_con_LSF_1_LE_p8(:,i));
    std_beta_con_LSF_1_LE_p8(1,i) = nanstd(Winstore_beta_con_LSF_1_LE_p8(:,i));
    std_dist_con_LSF_1_LE_p8(1,i) = nanstd(trackingerror_con_LSF_1_LE_p8(:,i));
    
    mean_x_con_LSF_1_LE_1(1,i) = nanmean(Winstore_x_con_LSF_1_LE_1(:,i));
    mean_y_con_LSF_1_LE_1(1,i) = nanmean(Winstore_y_con_LSF_1_LE_1(:,i));
    mean_theta_con_LSF_1_LE_1(1,i) = nanmean(Winstore_theta_con_LSF_1_LE_1(:,i));
    mean_phi_con_LSF_1_LE_1(1,i) = nanmean(Winstore_phi_con_LSF_1_LE_1(:,i));
    mean_xdot_con_LSF_1_LE_1(1,i) = nanmean(Winstore_xdot_con_LSF_1_LE_1(:,i));
    mean_ydot_con_LSF_1_LE_1(1,i) = nanmean(Winstore_ydot_con_LSF_1_LE_1(:,i));
    mean_thetadot_con_LSF_1_LE_1(1,i) = nanmean(Winstore_thetadot_con_LSF_1_LE_1(:,i));
    mean_phidot_con_LSF_1_LE_1(1,i) = nanmean(Winstore_phidot_con_LSF_1_LE_1(:,i));
    mean_beta_con_LSF_1_LE_1(1,i) = nanmean(Winstore_beta_con_LSF_1_LE_1(:,i));
    mean_dist_con_LSF_1_LE_1(1,i) = nanmean(trackingerror_con_LSF_1_LE_1(:,i));
    std_x_con_LSF_1_LE_1(1,i) = nanstd(Winstore_x_con_LSF_1_LE_1(:,i));
    std_y_con_LSF_1_LE_1(1,i) = nanstd(Winstore_y_con_LSF_1_LE_1(:,i));
    std_theta_con_LSF_1_LE_1(1,i) = nanstd(Winstore_theta_con_LSF_1_LE_1(:,i));
    std_phi_con_LSF_1_LE_1(1,i) = nanstd(Winstore_phi_con_LSF_1_LE_1(:,i));
    std_xdot_con_LSF_1_LE_1(1,i) = nanstd(Winstore_xdot_con_LSF_1_LE_1(:,i));
    std_ydot_con_LSF_1_LE_1(1,i) = nanstd(Winstore_ydot_con_LSF_1_LE_1(:,i));
    std_thetadot_con_LSF_1_LE_1(1,i) = nanstd(Winstore_thetadot_con_LSF_1_LE_1(:,i));
    std_phidot_con_LSF_1_LE_1(1,i) = nanstd(Winstore_phidot_con_LSF_1_LE_1(:,i));
    std_beta_con_LSF_1_LE_1(1,i) = nanstd(Winstore_beta_con_LSF_1_LE_1(:,i));
    std_dist_con_LSF_1_LE_1(1,i) = nanstd(trackingerror_con_LSF_1_LE_1(:,i));
    
    mean_x_con_LSF_1_LE_2(1,i) = nanmean(Winstore_x_con_LSF_1_LE_2(:,i));
    mean_y_con_LSF_1_LE_2(1,i) = nanmean(Winstore_y_con_LSF_1_LE_2(:,i));
    mean_theta_con_LSF_1_LE_2(1,i) = nanmean(Winstore_theta_con_LSF_1_LE_2(:,i));
    mean_phi_con_LSF_1_LE_2(1,i) = nanmean(Winstore_phi_con_LSF_1_LE_2(:,i));
    mean_xdot_con_LSF_1_LE_2(1,i) = nanmean(Winstore_xdot_con_LSF_1_LE_2(:,i));
    mean_ydot_con_LSF_1_LE_2(1,i) = nanmean(Winstore_ydot_con_LSF_1_LE_2(:,i));
    mean_thetadot_con_LSF_1_LE_2(1,i) = nanmean(Winstore_thetadot_con_LSF_1_LE_2(:,i));
    mean_phidot_con_LSF_1_LE_2(1,i) = nanmean(Winstore_phidot_con_LSF_1_LE_2(:,i));
    mean_beta_con_LSF_1_LE_2(1,i) = nanmean(Winstore_beta_con_LSF_1_LE_2(:,i));
    mean_dist_con_LSF_1_LE_2(1,i) = nanmean(trackingerror_con_LSF_1_LE_2(:,i));
    std_x_con_LSF_1_LE_2(1,i) = nanstd(Winstore_x_con_LSF_1_LE_2(:,i));
    std_y_con_LSF_1_LE_2(1,i) = nanstd(Winstore_y_con_LSF_1_LE_2(:,i));
    std_theta_con_LSF_1_LE_2(1,i) = nanstd(Winstore_theta_con_LSF_1_LE_2(:,i));
    std_phi_con_LSF_1_LE_2(1,i) = nanstd(Winstore_phi_con_LSF_1_LE_2(:,i));
    std_xdot_con_LSF_1_LE_2(1,i) = nanstd(Winstore_xdot_con_LSF_1_LE_2(:,i));
    std_ydot_con_LSF_1_LE_2(1,i) = nanstd(Winstore_ydot_con_LSF_1_LE_2(:,i));
    std_thetadot_con_LSF_1_LE_2(1,i) = nanstd(Winstore_thetadot_con_LSF_1_LE_2(:,i));
    std_phidot_con_LSF_1_LE_2(1,i) = nanstd(Winstore_phidot_con_LSF_1_LE_2(:,i));
    std_beta_con_LSF_1_LE_2(1,i) = nanstd(Winstore_beta_con_LSF_1_LE_2(:,i));
    std_dist_con_LSF_1_LE_2(1,i) = nanstd(trackingerror_con_LSF_1_LE_2(:,i));
    
    %numel(A(~isnan(A)))
    numel_con_LSF_1_LE_0(1,i) = numel(Winstore_x_con_LSF_1_LE_0(~isnan(Winstore_x_con_LSF_1_LE_0(:,i)))); 
    numel_con_LSF_1_LE_p2(1,i) = numel(Winstore_x_con_LSF_1_LE_p2(~isnan(Winstore_x_con_LSF_1_LE_p2(:,i)))); 
    numel_con_LSF_1_LE_p4(1,i) = numel(Winstore_x_con_LSF_1_LE_p4(~isnan(Winstore_x_con_LSF_1_LE_p4(:,i)))); 
    numel_con_LSF_1_LE_p6(1,i) = numel(Winstore_x_con_LSF_1_LE_p6(~isnan(Winstore_x_con_LSF_1_LE_p6(:,i)))); 
    numel_con_LSF_1_LE_p8(1,i) = numel(Winstore_x_con_LSF_1_LE_p8(~isnan(Winstore_x_con_LSF_1_LE_p8(:,i)))); 
    numel_con_LSF_1_LE_1(1,i) = numel(Winstore_x_con_LSF_1_LE_1(~isnan(Winstore_x_con_LSF_1_LE_1(:,i)))); 
    numel_con_LSF_1_LE_2(1,i) = numel(Winstore_x_con_LSF_1_LE_2(~isnan(Winstore_x_con_LSF_1_LE_2(:,i)))); 
end

for i = 1:hws
    mean_impcost_con_LSF_1_LE_0(1,i) = nanmean(cost_import_con_LSF_1_LE_0(:,i));
    mean_impcost_con_LSF_1_LE_p2(1,i) = nanmean(cost_import_con_LSF_1_LE_p2(:,i));
    mean_impcost_con_LSF_1_LE_p4(1,i) = nanmean(cost_import_con_LSF_1_LE_p4(:,i));
    mean_impcost_con_LSF_1_LE_p6(1,i) = nanmean(cost_import_con_LSF_1_LE_p6(:,i));
    mean_impcost_con_LSF_1_LE_p8(1,i) = nanmean(cost_import_con_LSF_1_LE_p8(:,i));
    mean_impcost_con_LSF_1_LE_1(1,i) = nanmean(cost_import_con_LSF_1_LE_1(:,i));
    mean_impcost_con_LSF_1_LE_2(1,i) = nanmean(cost_import_con_LSF_1_LE_2(:,i));
    
    std_impcost_con_LSF_1_LE_0(1,i) = nanstd(cost_import_con_LSF_1_LE_0(:,i));
    std_impcost_con_LSF_1_LE_p2(1,i) = nanstd(cost_import_con_LSF_1_LE_p2(:,i));
    std_impcost_con_LSF_1_LE_p4(1,i) = nanstd(cost_import_con_LSF_1_LE_p4(:,i));
    std_impcost_con_LSF_1_LE_p6(1,i) = nanstd(cost_import_con_LSF_1_LE_p6(:,i));
    std_impcost_con_LSF_1_LE_p8(1,i) = nanstd(cost_import_con_LSF_1_LE_p8(:,i));
    std_impcost_con_LSF_1_LE_1(1,i) = nanstd(cost_import_con_LSF_1_LE_1(:,i));
    std_impcost_con_LSF_1_LE_2(1,i) = nanstd(cost_import_con_LSF_1_LE_2(:,i));
end
toc

%% Extract F, alpha, and tau values
disp('Extracting F, alpha, and tau values') 
%Magnitude of forces
F_con_LSF_1_LE_0 = ICsave_am_con_LSF_1_LE_0(:,9:12:end);
F_con_LSF_1_LE_p2 = ICsave_am_con_LSF_1_LE_p2(:,9:12:end);
F_con_LSF_1_LE_p4 = ICsave_am_con_LSF_1_LE_p4(:,9:12:end);
F_con_LSF_1_LE_p6 = ICsave_am_con_LSF_1_LE_p6(:,9:12:end);
F_con_LSF_1_LE_p8 = ICsave_am_con_LSF_1_LE_p8(:,9:12:end);
F_con_LSF_1_LE_1 = ICsave_am_con_LSF_1_LE_1(:,9:12:end);
F_con_LSF_1_LE_2 = ICsave_am_con_LSF_1_LE_2(:,9:12:end);

%Direction of aforementioned forces
alpha_con_LSF_1_LE_0 = ICsave_am_con_LSF_1_LE_0(:,10:12:end);
alpha_con_LSF_1_LE_p2 = ICsave_am_con_LSF_1_LE_p2(:,10:12:end);
alpha_con_LSF_1_LE_p4 = ICsave_am_con_LSF_1_LE_p4(:,10:12:end);
alpha_con_LSF_1_LE_p6 = ICsave_am_con_LSF_1_LE_p6(:,10:12:end);
alpha_con_LSF_1_LE_p8 = ICsave_am_con_LSF_1_LE_p8(:,10:12:end);
alpha_con_LSF_1_LE_1 = ICsave_am_con_LSF_1_LE_1(:,10:12:end);
alpha_con_LSF_1_LE_2 = ICsave_am_con_LSF_1_LE_2(:,10:12:end);

%Torque exerted about the pin joint
tau_con_LSF_1_LE_0 = ICsave_am_con_LSF_1_LE_0(:,11:12:end);
tau_con_LSF_1_LE_p2 = ICsave_am_con_LSF_1_LE_p2(:,11:12:end);
tau_con_LSF_1_LE_p4 = ICsave_am_con_LSF_1_LE_p4(:,11:12:end);
tau_con_LSF_1_LE_p6 = ICsave_am_con_LSF_1_LE_p6(:,11:12:end);
tau_con_LSF_1_LE_p8 = ICsave_am_con_LSF_1_LE_p8(:,11:12:end);
tau_con_LSF_1_LE_1 = ICsave_am_con_LSF_1_LE_1(:,11:12:end);
tau_con_LSF_1_LE_2 = ICsave_am_con_LSF_1_LE_2(:,11:12:end);

%Remove the pesky zeros
F_con_LSF_1_LE_0(F_con_LSF_1_LE_0==0) = NaN;
F_con_LSF_1_LE_p2(F_con_LSF_1_LE_p2==0) = NaN;
F_con_LSF_1_LE_p4(F_con_LSF_1_LE_p4==0) = NaN;
F_con_LSF_1_LE_p6(F_con_LSF_1_LE_p6==0) = NaN;
F_con_LSF_1_LE_p8(F_con_LSF_1_LE_p8==0) = NaN;
F_con_LSF_1_LE_1(F_con_LSF_1_LE_1==0) = NaN;
F_con_LSF_1_LE_2(F_con_LSF_1_LE_2==0) = NaN;

alpha_con_LSF_1_LE_0(alpha_con_LSF_1_LE_0==0) = NaN;
alpha_con_LSF_1_LE_p2(alpha_con_LSF_1_LE_p2==0) = NaN;
alpha_con_LSF_1_LE_p4(alpha_con_LSF_1_LE_p4==0) = NaN;
alpha_con_LSF_1_LE_p6(alpha_con_LSF_1_LE_p6==0) = NaN;
alpha_con_LSF_1_LE_p8(alpha_con_LSF_1_LE_p8==0) = NaN;
alpha_con_LSF_1_LE_1(alpha_con_LSF_1_LE_1==0) = NaN;
alpha_con_LSF_1_LE_2(alpha_con_LSF_1_LE_2==0) = NaN;

tau_con_LSF_1_LE_0(tau_con_LSF_1_LE_0==0) = NaN;
tau_con_LSF_1_LE_p2(tau_con_LSF_1_LE_p2==0) = NaN;
tau_con_LSF_1_LE_p4(tau_con_LSF_1_LE_p4==0) = NaN;
tau_con_LSF_1_LE_p6(tau_con_LSF_1_LE_p6==0) = NaN;
tau_con_LSF_1_LE_p8(tau_con_LSF_1_LE_p8==0) = NaN;
tau_con_LSF_1_LE_1(tau_con_LSF_1_LE_1==0) = NaN;
tau_con_LSF_1_LE_2(tau_con_LSF_1_LE_2==0) = NaN;

disp('Done curating, now save the files')

%% Save the files

%con - LSF_1_LE_0
save('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_0.mat','mean_x_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_0.mat','mean_y_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_0.mat','mean_theta_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_0.mat','mean_phi_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_0.mat','mean_xdot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_0.mat','mean_ydot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_0.mat','mean_thetadot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_0.mat','mean_phidot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_0.mat','mean_beta_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_0.mat','mean_dist_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_0.mat','std_x_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_0.mat','std_y_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_0.mat','std_theta_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_0.mat','std_phi_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_0.mat','std_xdot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_0.mat','std_ydot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_0.mat','std_thetadot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_0.mat','std_phidot_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_0.mat','std_beta_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_0.mat','std_dist_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_0.mat','F_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_0.mat','alpha_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_0.mat','tau_con_LSF_1_LE_0');

%con - LSF_1_LE_p2
save('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p2.mat','mean_x_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p2.mat','mean_y_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p2.mat','mean_theta_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p2.mat','mean_phi_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p2.mat','mean_xdot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p2.mat','mean_ydot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p2.mat','mean_thetadot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p2.mat','mean_phidot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p2.mat','mean_beta_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p2.mat','mean_dist_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p2.mat','std_x_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p2.mat','std_y_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p2.mat','std_theta_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p2.mat','std_phi_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p2.mat','std_xdot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p2.mat','std_ydot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p2.mat','std_thetadot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p2.mat','std_phidot_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p2.mat','std_beta_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p2.mat','std_dist_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p2.mat','F_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p2.mat','alpha_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p2.mat','tau_con_LSF_1_LE_p2');

%con - LSF_1_LE_p4
save('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p4.mat','mean_x_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p4.mat','mean_y_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p4.mat','mean_theta_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p4.mat','mean_phi_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p4.mat','mean_xdot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p4.mat','mean_ydot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p4.mat','mean_thetadot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p4.mat','mean_phidot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p4.mat','mean_beta_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p4.mat','mean_dist_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p4.mat','std_x_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p4.mat','std_y_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p4.mat','std_theta_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p4.mat','std_phi_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p4.mat','std_xdot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p4.mat','std_ydot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p4.mat','std_thetadot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p4.mat','std_phidot_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p4.mat','std_beta_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p4.mat','std_dist_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p4.mat','F_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p4.mat','alpha_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p4.mat','tau_con_LSF_1_LE_p4');

%con - LSF_1_LE_p6
save('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p6.mat','mean_x_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p6.mat','mean_y_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p6.mat','mean_theta_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p6.mat','mean_phi_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p6.mat','mean_xdot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p6.mat','mean_ydot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p6.mat','mean_thetadot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p6.mat','mean_phidot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p6.mat','mean_beta_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p6.mat','mean_dist_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p6.mat','std_x_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p6.mat','std_y_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p6.mat','std_theta_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p6.mat','std_phi_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p6.mat','std_xdot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p6.mat','std_ydot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p6.mat','std_thetadot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p6.mat','std_phidot_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p6.mat','std_beta_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p6.mat','std_dist_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p6.mat','F_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p6.mat','alpha_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p6.mat','tau_con_LSF_1_LE_p6');

%con - LSF_1_LE_p8
save('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_p8.mat','mean_x_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_p8.mat','mean_y_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_p8.mat','mean_theta_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_p8.mat','mean_phi_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_p8.mat','mean_xdot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_p8.mat','mean_ydot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_p8.mat','mean_thetadot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_p8.mat','mean_phidot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_p8.mat','mean_beta_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_p8.mat','mean_dist_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_p8.mat','std_x_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_p8.mat','std_y_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_p8.mat','std_theta_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_p8.mat','std_phi_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_p8.mat','std_xdot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_p8.mat','std_ydot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_p8.mat','std_thetadot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_p8.mat','std_phidot_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_p8.mat','std_beta_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_p8.mat','std_dist_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p8.mat','F_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p8.mat','alpha_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p8.mat','tau_con_LSF_1_LE_p8');

%con - LSF_1_LE_1
save('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_1.mat','mean_x_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_1.mat','mean_y_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_1.mat','mean_theta_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_1.mat','mean_phi_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_1.mat','mean_xdot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_1.mat','mean_ydot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_1.mat','mean_thetadot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_1.mat','mean_phidot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_1.mat','mean_beta_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_1.mat','mean_dist_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_1.mat','std_x_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_1.mat','std_y_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_1.mat','std_theta_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_1.mat','std_phi_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_1.mat','std_xdot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_1.mat','std_ydot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_1.mat','std_thetadot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_1.mat','std_phidot_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_1.mat','std_beta_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_1.mat','std_dist_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_1.mat','F_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_1.mat','alpha_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_1.mat','tau_con_LSF_1_LE_1');

%con - LSF_1_LE_2
save('../CuratedData_MPC/LSF_1/mean_x_con_LSF_1_LE_2.mat','mean_x_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_y_con_LSF_1_LE_2.mat','mean_y_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_theta_con_LSF_1_LE_2.mat','mean_theta_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_phi_con_LSF_1_LE_2.mat','mean_phi_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_xdot_con_LSF_1_LE_2.mat','mean_xdot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_ydot_con_LSF_1_LE_2.mat','mean_ydot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_thetadot_con_LSF_1_LE_2.mat','mean_thetadot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_phidot_con_LSF_1_LE_2.mat','mean_phidot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_beta_con_LSF_1_LE_2.mat','mean_beta_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/mean_dist_con_LSF_1_LE_2.mat','mean_dist_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_x_con_LSF_1_LE_2.mat','std_x_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_y_con_LSF_1_LE_2.mat','std_y_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_theta_con_LSF_1_LE_2.mat','std_theta_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_phi_con_LSF_1_LE_2.mat','std_phi_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_xdot_con_LSF_1_LE_2.mat','std_xdot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_ydot_con_LSF_1_LE_2.mat','std_ydot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_thetadot_con_LSF_1_LE_2.mat','std_thetadot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_phidot_con_LSF_1_LE_2.mat','std_phidot_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_beta_con_LSF_1_LE_2.mat','std_beta_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/std_dist_con_LSF_1_LE_2.mat','std_dist_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_2.mat','F_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_2.mat','alpha_con_LSF_1_LE_2');
save('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_2.mat','tau_con_LSF_1_LE_2');

%numel
save('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_0.mat','numel_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p2.mat','numel_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p4.mat','numel_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p6.mat','numel_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_p8.mat','numel_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_1.mat','numel_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/numel_con_LSF_1_LE_2.mat','numel_con_LSF_1_LE_2');

%imported cost
save('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_0.mat','mean_impcost_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p2.mat','mean_impcost_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p4.mat','mean_impcost_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p6.mat','mean_impcost_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_p8.mat','mean_impcost_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_1.mat','mean_impcost_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/mean_impcost_con_LSF_1_LE_2.mat','mean_impcost_con_LSF_1_LE_2');

save('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_0.mat','std_impcost_con_LSF_1_LE_0');
save('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p2.mat','std_impcost_con_LSF_1_LE_p2');
save('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p4.mat','std_impcost_con_LSF_1_LE_p4');
save('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p6.mat','std_impcost_con_LSF_1_LE_p6');
save('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_p8.mat','std_impcost_con_LSF_1_LE_p8');
save('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_1.mat','std_impcost_con_LSF_1_LE_1');
save('../CuratedData_MPC/LSF_1/std_impcost_con_LSF_1_LE_2.mat','std_impcost_con_LSF_1_LE_2');


%tepfr - tracking error per full run
save('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_0','mean_tepfr_con_LSF_1_LE_0')
save('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p2','mean_tepfr_con_LSF_1_LE_p2')
save('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p4','mean_tepfr_con_LSF_1_LE_p4')
save('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p6','mean_tepfr_con_LSF_1_LE_p6')
save('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_p8','mean_tepfr_con_LSF_1_LE_p8')
save('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_1','mean_tepfr_con_LSF_1_LE_1')
save('../CuratedData_MPC/LSF_1/mean_tepfr_con_LSF_1_LE_2','mean_tepfr_con_LSF_1_LE_2')

%Mean of cost per full run -- for stats purposes
save('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_0','mean_cost_con_LSF_1_LE_0')
save('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p2','mean_cost_con_LSF_1_LE_p2')
save('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p4','mean_cost_con_LSF_1_LE_p4')
save('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p6','mean_cost_con_LSF_1_LE_p6')
save('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_p8','mean_cost_con_LSF_1_LE_p8')
save('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_1','mean_cost_con_LSF_1_LE_1')
save('../CuratedData_MPC/LSF_1/mean_cost_con_LSF_1_LE_2','mean_cost_con_LSF_1_LE_2')

%% Done, son
disp('Done saving the files')
disp(['Time at end: ',datestr(now)])