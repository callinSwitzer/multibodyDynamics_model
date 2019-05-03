%7/9/2018 Script developed to plot the F, alpha, tau distributions
%8/23/2018 Script modified for LengthScaleFactor (LSF)
%10/23/2018 Script modified for size extension (LE)
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

%% List the directory stuff
listOfconLSF_1_LE_0_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_0.mat'); 
listOfconLSF_1_LE_p2_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p2.mat'); 
listOfconLSF_1_LE_p4_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p4.mat'); 
listOfconLSF_1_LE_p6_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p6.mat'); 
listOfconLSF_1_LE_p8_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_p8.mat'); 
listOfconLSF_1_LE_1_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_1.mat'); 
listOfconLSF_1_LE_2_Files = dir('../SimData_MPC/LSF_1/Winstore_*_sp_con_LSF_1_LE_2.mat'); 

LSF_1 = 1;

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
PartPath = 2000;
timestep = numel(Tstore(1:(end-100)))/hws;
tderiv = Tstore(2);
t_end = Tstore(end-100);

% m = (rho*(4/3)*pi*(bhead^2)*ahead + rho*(4/3)*pi*(bbutt^2)*abutt);

%Create the sum of primes signal
signal_amp = 5; %in cm 
prime_f = [0.2, 0.3, 0.5, 0.7, 1.1, 1.7, 2.9, 4.3, 7.9, 13.7, 19.9]; %in Hz
prime_a = (signal_amp./(2*pi.*prime_f)).*(2*pi.*prime_f(1)); %in cm
prime_ph = zeros(1,numel(prime_f)); %Not sure if we'll require a phase

%Goal criteria (for the cost function)
y_g = prime_a(1)*sin(2*pi*prime_f(1)*Tstore(1:(end-100)) + prime_ph(1)) +...
    prime_a(2)*sin(2*pi*prime_f(2)*Tstore(1:(end-100)) + prime_ph(2)) +...
    prime_a(3)*sin(2*pi*prime_f(3)*Tstore(1:(end-100)) + prime_ph(3)) +...
    prime_a(4)*sin(2*pi*prime_f(4)*Tstore(1:(end-100)) + prime_ph(4)) +...
    prime_a(5)*sin(2*pi*prime_f(5)*Tstore(1:(end-100)) + prime_ph(5)) +...
    prime_a(6)*sin(2*pi*prime_f(6)*Tstore(1:(end-100)) + prime_ph(6)) +...
    prime_a(7)*sin(2*pi*prime_f(7)*Tstore(1:(end-100)) + prime_ph(7)) +...
    prime_a(8)*sin(2*pi*prime_f(8)*Tstore(1:(end-100)) + prime_ph(8)) +...
    prime_a(9)*sin(2*pi*prime_f(9)*Tstore(1:(end-100)) + prime_ph(9)) +...
    prime_a(10)*sin(2*pi*prime_f(10)*Tstore(1:(end-100)) + prime_ph(10)) +...
    prime_a(11)*sin(2*pi*prime_f(11)*Tstore(1:(end-100)) + prime_ph(11));

ydot_g = 2*pi*prime_a(1)*prime_f(1)*cos(2*pi*prime_f(1)*Tstore(1:(end-100)) + prime_ph(1)) +...
    2*pi*prime_a(2)*prime_f(2)*cos(2*pi*prime_f(2)*Tstore(1:(end-100)) + prime_ph(2)) +...
    2*pi*prime_a(3)*prime_f(3)*cos(2*pi*prime_f(3)*Tstore(1:(end-100)) + prime_ph(3)) +...
    2*pi*prime_a(4)*prime_f(4)*cos(2*pi*prime_f(4)*Tstore(1:(end-100)) + prime_ph(4)) +...
    2*pi*prime_a(5)*prime_f(5)*cos(2*pi*prime_f(5)*Tstore(1:(end-100)) + prime_ph(5)) +...
    2*pi*prime_a(6)*prime_f(6)*cos(2*pi*prime_f(6)*Tstore(1:(end-100)) + prime_ph(6)) +...
    2*pi*prime_a(7)*prime_f(7)*cos(2*pi*prime_f(7)*Tstore(1:(end-100)) + prime_ph(7)) +...
    2*pi*prime_a(8)*prime_f(8)*cos(2*pi*prime_f(8)*Tstore(1:(end-100)) + prime_ph(8)) +...
    2*pi*prime_a(9)*prime_f(9)*cos(2*pi*prime_f(9)*Tstore(1:(end-100)) + prime_ph(9)) +...
    2*pi*prime_a(10)*prime_f(10)*cos(2*pi*prime_f(10)*Tstore(1:(end-100)) + prime_ph(10)) +...
    2*pi*prime_a(11)*prime_f(11)*cos(2*pi*prime_f(11)*Tstore(1:(end-100)) + prime_ph(11));

y_g = y_g';
x_g = zeros(1,numel(y_g));
theta_g = zeros(1,numel(y_g));
theta_g(1,:) = pi/4;

ydot_g = ydot_g';
xdot_g = zeros(1,numel(y_g));
thetadot_g = zeros(1,numel(y_g));

time_hws(1:PartPath) = Tstore(1:(timestep*hws/PartPath):(end-100));

NumOfconLSF_1_LE_0_Files = numel(listOfconLSF_1_LE_0_Files);
NumOfconLSF_1_LE_p2_Files = numel(listOfconLSF_1_LE_p2_Files);
NumOfconLSF_1_LE_p4_Files = numel(listOfconLSF_1_LE_p4_Files);
NumOfconLSF_1_LE_p6_Files = numel(listOfconLSF_1_LE_p6_Files);
NumOfconLSF_1_LE_p8_Files = numel(listOfconLSF_1_LE_p8_Files);
NumOfconLSF_1_LE_1_Files = numel(listOfconLSF_1_LE_1_Files);
NumOfconLSF_1_LE_2_Files = numel(listOfconLSF_1_LE_2_Files);

%% Import all the variables
disp('Importing the relevant data')

%Magnitude of applied forces
load('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/F_con_LSF_1_LE_2.mat')

%Direction of applied forces
load('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/alpha_con_LSF_1_LE_2.mat')

%Torque applied about the pin joint
load('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_0.mat')
load('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p6.mat')
load('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_p8.mat')
load('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_1.mat')
load('../CuratedData_MPC/LSF_1/tau_con_LSF_1_LE_2.mat')

%% Find the min & max values for...
%Forces
Fmax_vec = [max(max(F_con_LSF_1_LE_0)),...
    max(max(F_con_LSF_1_LE_p2)),...
    max(max(F_con_LSF_1_LE_p4)),...
    max(max(F_con_LSF_1_LE_p6)),...
    max(max(F_con_LSF_1_LE_p8)),...
    max(max(F_con_LSF_1_LE_1)),...
    max(max(F_con_LSF_1_LE_2))];
F_max = max(Fmax_vec);

Fmin_vec = [min(min(F_con_LSF_1_LE_0)),...
    min(min(F_con_LSF_1_LE_p2)),...
    min(min(F_con_LSF_1_LE_p4)),...
    min(min(F_con_LSF_1_LE_p6)),...
    min(min(F_con_LSF_1_LE_p8)),...
    min(min(F_con_LSF_1_LE_1)),...
    min(min(F_con_LSF_1_LE_2))];
F_min = min(Fmin_vec);

%Tau
tauMax_vec = [max(max(tau_con_LSF_1_LE_0)),...
    max(max(tau_con_LSF_1_LE_p2)),...
    max(max(tau_con_LSF_1_LE_p4)),...
    max(max(tau_con_LSF_1_LE_p6)),...
    max(max(tau_con_LSF_1_LE_p8)),...
    max(max(tau_con_LSF_1_LE_1)),...
    max(max(tau_con_LSF_1_LE_2))];
tau_max = max(tauMax_vec);

tauMin_vec = [min(min(tau_con_LSF_1_LE_0)),...
    min(min(tau_con_LSF_1_LE_p2)),...
    min(min(tau_con_LSF_1_LE_p4)),...
    min(min(tau_con_LSF_1_LE_p6)),...
    min(min(tau_con_LSF_1_LE_p8)),...
    min(min(tau_con_LSF_1_LE_1)),...
    min(min(tau_con_LSF_1_LE_2))];
tau_min = min(tauMin_vec);

%% Magnitude of applied force plot
fig1 = figure(1);

%Note: units of force are in g*cm/(s^2) OR 10^(-5) Newtons
%y-motion
subplot(8,1,1)
plot(Tstore(1:(end-100)),y_g,'-k','LineWidth',2)
title(['LSF: ',num2str(LSF_1)]);
ylabel('y goal (cm)'); 
xlabel('time (s)');
axis([0, Tstore(end), -11, 11])

%con - LSF_1_LE_0
subplot(8,1,2)
plot(time_hws,F_con_LSF_1_LE_0./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['LE: 0, (n = ',num2str(NumOfconLSF_1_LE_0_Files),')']);
axis([0, Tstore(end), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p2
subplot(8,1,3)
plot(time_hws,F_con_LSF_1_LE_p2./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['LE: 0.2, (n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']);
axis([0, Tstore(end), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p4
subplot(8,1,4)
plot(time_hws,F_con_LSF_1_LE_p4./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['LE: 0.4, (n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']);
axis([0, Tstore(end), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p6
subplot(8,1,5)
plot(time_hws,F_con_LSF_1_LE_p6./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['LE: 0.6, (n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']);
axis([0, Tstore(end), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p8
subplot(8,1,6)
plot(time_hws,F_con_LSF_1_LE_p8./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['LE: 0.8, (n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']);
axis([0, Tstore(end), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_1
subplot(8,1,7)
plot(time_hws,F_con_LSF_1_LE_1./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['LE: 1, (n = ',num2str(NumOfconLSF_1_LE_1_Files),')']);
axis([0, Tstore(end), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_2
subplot(8,1,8)
plot(time_hws,F_con_LSF_1_LE_2./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['LE: 2, (n = ',num2str(NumOfconLSF_1_LE_2_Files),')']);
axis([0, Tstore(end), F_min./(1000*100), F_max./(1000*100)])

set(fig1, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/LSF_1/Fig9_F_LSF_1_LE_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig9_F_LSF_1_LE_v12h.jpg');

disp('Figure 1 done')

%% Direction of applied force plot
fig2 = figure(2);

%y-motion
subplot(8,1,1)
plot(Tstore(1:(end-100)),y_g,'-k','LineWidth',2)
title(['LSF: ',num2str(LSF_1)]);
ylabel('y goal (cm)');  
xlabel('time (s)');
axis([0, Tstore(end), -11, 11])

%con - LSF_1_LE_0
subplot(8,1,2)
plot(time_hws,alpha_con_LSF_1_LE_0,'ok','MarkerSize',2)
ylabel('Alpha (rad)'); 
xlabel(['LE: 0, (n = ',num2str(NumOfconLSF_1_LE_0_Files),')']);
axis([0, Tstore(end), 0, 2*pi])

%con - LSF_1_LE_p2
subplot(8,1,3)
plot(time_hws,alpha_con_LSF_1_LE_p2,'ok','MarkerSize',2)
ylabel('Alpha (rad)'); 
xlabel(['LE: 0.2, (n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']);
axis([0, Tstore(end), 0, 2*pi])

%con - LSF_1_LE_p4
subplot(8,1,4)
plot(time_hws,alpha_con_LSF_1_LE_p4,'ok','MarkerSize',2)
ylabel('Alpha (rad)'); 
xlabel(['LE: 0.4, (n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']);
axis([0, Tstore(end), 0, 2*pi])

%con - LSF_1_LE_p6
subplot(8,1,5)
plot(time_hws,alpha_con_LSF_1_LE_p6,'ok','MarkerSize',2)
ylabel('Alpha (rad)'); 
xlabel(['LE: 0.6, (n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']);
axis([0, Tstore(end), 0, 2*pi])

%con - LSF_1_LE_p8
subplot(8,1,6)
plot(time_hws,alpha_con_LSF_1_LE_p8,'ok','MarkerSize',2)
ylabel('Alpha (rad)'); 
xlabel(['LE: 0.8, (n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']);
axis([0, Tstore(end), 0, 2*pi])

%con - LSF_1_LE_1
subplot(8,1,7)
plot(time_hws,alpha_con_LSF_1_LE_1,'ok','MarkerSize',2)
ylabel('Alpha (rad)'); 
xlabel(['LE: 1, (n = ',num2str(NumOfconLSF_1_LE_1_Files),')']);
axis([0, Tstore(end), 0, 2*pi])

%con - LSF_1_LE_2
subplot(8,1,8)
plot(time_hws,alpha_con_LSF_1_LE_2,'ok','MarkerSize',2)
ylabel('Alpha (rad)'); 
xlabel(['LE: 2, (n = ',num2str(NumOfconLSF_1_LE_2_Files),')']);
axis([0, Tstore(end), 0, 2*pi])

set(fig2, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/LSF_1/Fig10_alpha_LSF_1_LE_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig10_alpha_LSF_1_LE_v12h.jpg');

disp('Figure 2 done')

%% Direction of applied torque plot
fig3 = figure(3);

%Note: units of torque are in g*(cm^2)/(s^2) 

zeroline(1:numel(time_hws)) = 0;

%y-motion
subplot(8,1,1)
plot(Tstore(1:(end-100)),y_g,'-k','LineWidth',2)
ylabel('y goal (cm)'); 
title(['LSF: ',num2str(LSF_1)]);
xlabel('time (s)');
axis([0, Tstore(end), -11, 11])

%con - LSF_1_LE_0
subplot(8,1,2)
plot(time_hws,tau_con_LSF_1_LE_0./(100^2*1000),'ok','MarkerSize',2)
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel('Tau (N*m)'); 
xlabel(['LE: 0, (n = ',num2str(NumOfconLSF_1_LE_0_Files),')']);
axis([0, Tstore(end), tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p2
subplot(8,1,3)
plot(time_hws,tau_con_LSF_1_LE_p2./(100^2*1000),'ok','MarkerSize',2)
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel('Tau (N*m)'); 
xlabel(['LE: 0.2, (n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']);
axis([0, Tstore(end), tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p4
subplot(8,1,4)
plot(time_hws,tau_con_LSF_1_LE_p4./(100^2*1000),'ok','MarkerSize',2)
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel('Tau (N*m)'); 
xlabel(['LE: 0.4, (n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']);
axis([0, Tstore(end), tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p6
subplot(8,1,5)
plot(time_hws,tau_con_LSF_1_LE_p6./(100^2*1000),'ok','MarkerSize',2)
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel('Tau (N*m)'); 
xlabel(['LE: 0.6, (n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']);
axis([0, Tstore(end), tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p8
subplot(8,1,6)
plot(time_hws,tau_con_LSF_1_LE_p8./(100^2*1000),'ok','MarkerSize',2)
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel('Tau (N*m)'); 
xlabel(['LE: 0.8, (n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']);
axis([0, Tstore(end), tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_1
subplot(8,1,7)
plot(time_hws,tau_con_LSF_1_LE_1./(100^2*1000),'ok','MarkerSize',2)
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel('Tau (N*m)'); 
xlabel(['LE: 1, (n = ',num2str(NumOfconLSF_1_LE_1_Files),')']);
axis([0, Tstore(end), tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_2
subplot(8,1,8)
plot(time_hws,tau_con_LSF_1_LE_2./(100^2*1000),'ok','MarkerSize',2)
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel('Tau (N*m)'); 
xlabel(['LE: 2, (n = ',num2str(NumOfconLSF_1_LE_2_Files),')']);
axis([0, Tstore(end), tau_min./(100^2*1000), tau_max./(100^2*1000)])

set(fig3, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/LSF_1/Fig11_tau_LSF_1_LE_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig11_tau_LSF_1_LE_v12h.jpg');

disp('Figure 3 done')

%% Force w.r.t. alpha
fig4 = figure(4);

%Note: units of force are in g*cm/(s^2) OR 10^(-5) Newtons

%con - LSF_1_LE_0
subplot(7,1,1)
title(['LSF: ',num2str(LSF_1)]);
plot(alpha_con_LSF_1_LE_0,F_con_LSF_1_LE_0./(1000*100),'ok','MarkerSize',2)
title(['LSF: ',num2str(LSF_1)]);
ylabel('Force (N)'); 
xlabel(['Alpha (rad), LE: 0, (n = ',num2str(NumOfconLSF_1_LE_0_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p2
subplot(7,1,2)
plot(alpha_con_LSF_1_LE_p2,F_con_LSF_1_LE_p2./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Alpha (rad), LE: 0.2, (n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p4
subplot(7,1,3)
plot(alpha_con_LSF_1_LE_p4,F_con_LSF_1_LE_p4./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Alpha (rad), LE: 0.4, (n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p6
subplot(7,1,4)
plot(alpha_con_LSF_1_LE_p6,F_con_LSF_1_LE_p6./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Alpha (rad), LE: 0.6, (n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p8
subplot(7,1,5)
plot(alpha_con_LSF_1_LE_p8,F_con_LSF_1_LE_p8./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Alpha (rad), LE: 0.8, (n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_1
subplot(7,1,6)
plot(alpha_con_LSF_1_LE_1,F_con_LSF_1_LE_1./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Alpha (rad), LE: 1, (n = ',num2str(NumOfconLSF_1_LE_1_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_2
subplot(7,1,7)
plot(alpha_con_LSF_1_LE_2,F_con_LSF_1_LE_2./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Alpha (rad), LE: 2, (n = ',num2str(NumOfconLSF_1_LE_2_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100)])

set(fig4, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/LSF_1/Fig12_F_alpha_LSF_1_LE_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig12_F_alpha_LSF_1_LE_v12h.jpg');

disp('Figure 4 done')

%% Applied force w.r.t. applied torque
fig5 = figure(5);

%Note: units of force are in g*cm/(s^2) OR 10^(-5) Newtons
%con - LSF_1_LE_0
subplot(7,1,1)
plot(tau_con_LSF_1_LE_0./(100^2*1000),F_con_LSF_1_LE_0./(1000*100),'ok','MarkerSize',2)
title(['LSF: ',num2str(LSF_1)]);
ylabel('Force (N)'); 
xlabel(['Tau (N*m), LE: 0, (n = ',num2str(NumOfconLSF_1_LE_0_Files),')']);
axis([tau_min./(100^2*1000), tau_max./(100^2*1000), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p2
subplot(7,1,2)
plot(tau_con_LSF_1_LE_p2./(100^2*1000),F_con_LSF_1_LE_p2./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Tau (N*m), LE: 0.2, (n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']);
axis([tau_min./(100^2*1000), tau_max./(100^2*1000), F_min./(1000*100), F_max./(1000*100)])


%con - LSF_1_LE_p4
subplot(7,1,3)
plot(tau_con_LSF_1_LE_p4./(100^2*1000),F_con_LSF_1_LE_p4./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Tau (N*m), LE: 0.4, (n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']);
axis([tau_min./(100^2*1000), tau_max./(100^2*1000), F_min./(1000*100), F_max./(1000*100)])


%con - LSF_1_LE_p6
subplot(7,1,4)
plot(tau_con_LSF_1_LE_p6./(100^2*1000),F_con_LSF_1_LE_p6./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Tau (N*m), LE: 0.6, (n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']);
axis([tau_min./(100^2*1000), tau_max./(100^2*1000), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_p8
subplot(7,1,5)
plot(tau_con_LSF_1_LE_p8./(100^2*1000),F_con_LSF_1_LE_p8./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Tau (N*m), LE: 0.8, (n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']);
axis([tau_min./(100^2*1000), tau_max./(100^2*1000), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_1
subplot(7,1,6)
plot(tau_con_LSF_1_LE_1./(100^2*1000),F_con_LSF_1_LE_1./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Tau (N*m), LE: 1, (n = ',num2str(NumOfconLSF_1_LE_1_Files),')']);
axis([tau_min./(100^2*1000), tau_max./(100^2*1000), F_min./(1000*100), F_max./(1000*100)])

%con - LSF_1_LE_2
subplot(7,1,7)
plot(tau_con_LSF_1_LE_2./(100^2*1000),F_con_LSF_1_LE_2./(1000*100),'ok','MarkerSize',2)
ylabel('Force (N)'); 
xlabel(['Tau (N*m), LE: 2, (n = ',num2str(NumOfconLSF_1_LE_2_Files),')']);
axis([tau_min./(100^2*1000), tau_max./(100^2*1000), F_min./(1000*100), F_max./(1000*100)])

set(fig5, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/LSF_1/Fig13_F_tau_LSF_1_LE_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig13_F_tau_LSF_1_LE_v12h.jpg');

disp('Figure 5 done')

%% Tau w.r.t. alpha
fig6 = figure(6);

%Note: units of force are in g*cm/(s^2) OR 10^(-5) Newtons
%con - LSF_1_LE_0
subplot(7,1,1)
plot(alpha_con_LSF_1_LE_0,tau_con_LSF_1_LE_0./(100^2*1000),'ok','MarkerSize',2)
title(['LSF: ',num2str(LSF_1)]);
ylabel('Tau (N*m)'); 
xlabel(['Alpha (rad), LE: 0, (n = ',num2str(NumOfconLSF_1_LE_0_Files),')']);
axis([0, 2*pi, tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p2
subplot(7,1,2)
plot(alpha_con_LSF_1_LE_p2,tau_con_LSF_1_LE_p2./(100^2*1000),'ok','MarkerSize',2)
ylabel('Tau (N*m)'); 
xlabel(['Alpha (rad), LE: 0.2, (n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']);
axis([0, 2*pi, tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p4
subplot(7,1,3)
plot(alpha_con_LSF_1_LE_p4,tau_con_LSF_1_LE_p4./(100^2*1000),'ok','MarkerSize',2)
ylabel('Tau (N*m)'); 
xlabel(['Alpha (rad), LE: 0.4, (n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']);
axis([0, 2*pi, tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p6
subplot(7,1,4)
plot(alpha_con_LSF_1_LE_p6,tau_con_LSF_1_LE_p6./(100^2*1000),'ok','MarkerSize',2)
ylabel('Tau (N*m)'); 
xlabel(['Alpha (rad), LE: 0.6, (n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']);
axis([0, 2*pi, tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p8
subplot(7,1,5)
plot(alpha_con_LSF_1_LE_p8,tau_con_LSF_1_LE_p8./(100^2*1000),'ok','MarkerSize',2)
ylabel('Tau (N*m)'); 
xlabel(['Alpha (rad), LE: 0.8, (n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']);
axis([0, 2*pi, tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_1
subplot(7,1,6)
plot(alpha_con_LSF_1_LE_1,tau_con_LSF_1_LE_1./(100^2*1000),'ok','MarkerSize',2)
ylabel('Tau (N*m)'); 
xlabel(['Alpha (rad), LE: 1, (n = ',num2str(NumOfconLSF_1_LE_1_Files),')']);
axis([0, 2*pi, tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_2
subplot(7,1,7)
plot(alpha_con_LSF_1_LE_2,tau_con_LSF_1_LE_2./(100^2*1000),'ok','MarkerSize',2)
ylabel('Tau (N*m)'); 
xlabel(['Alpha (rad), LE: 2, (n = ',num2str(NumOfconLSF_1_LE_2_Files),')']);
axis([0, 2*pi, tau_min./(100^2*1000), tau_max./(100^2*1000)])

set(fig6, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/LSF_1/Fig14_tau_alpha_LSF_1_LE_v12h.fig');
saveas(gcf,'../Figures_MPC/LSF_1/Fig14_tau_alpha_LSF_1_LE_v12h.jpg');

disp('Figure 6 done')

%% 3 axis plot - Force, alpha, tau
fig7 = figure(7);

%Note: units of force are in g*cm/(s^2) OR 10^(-5) Newtons
%con - LSF_1_LE_0
subplot(3,3,1)
plot3(alpha_con_LSF_1_LE_0,F_con_LSF_1_LE_0./(1000*100),...
    tau_con_LSF_1_LE_0./(100^2*1000),'ok','MarkerSize',2)
zlabel('Tau (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), LE: 0, (n = ',num2str(NumOfconLSF_1_LE_0_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100),...
    tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p2
subplot(3,3,2)
plot3(alpha_con_LSF_1_LE_p2,F_con_LSF_1_LE_p2./(1000*100),...
    tau_con_LSF_1_LE_p2./(100^2*1000),'ok','MarkerSize',2)
title(['LSF: ',num2str(LSF_1)]);
zlabel('Tau (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), LE: 0.2, (n = ',num2str(NumOfconLSF_1_LE_p2_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100),...
    tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p4
subplot(3,3,3)
plot3(alpha_con_LSF_1_LE_p4,F_con_LSF_1_LE_p4./(1000*100),...
    tau_con_LSF_1_LE_p4./(100^2*1000),'ok','MarkerSize',2)
zlabel('Tau (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), LE: 0.4, (n = ',num2str(NumOfconLSF_1_LE_p4_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100),...
    tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p6
subplot(3,3,4)
plot3(alpha_con_LSF_1_LE_p6,F_con_LSF_1_LE_p6./(1000*100),...
    tau_con_LSF_1_LE_p6./(100^2*1000),'ok','MarkerSize',2)
zlabel('Tau (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), LE: 0.6, (n = ',num2str(NumOfconLSF_1_LE_p6_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100),...
    tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_p8
subplot(3,3,5)
plot3(alpha_con_LSF_1_LE_p8,F_con_LSF_1_LE_p8./(1000*100),...
    tau_con_LSF_1_LE_p8./(100^2*1000),'ok','MarkerSize',2)
zlabel('Tau (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), LE: 0.8, (n = ',num2str(NumOfconLSF_1_LE_p8_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100),...
    tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_1
subplot(3,3,6)
plot3(alpha_con_LSF_1_LE_1,F_con_LSF_1_LE_1./(1000*100),...
    tau_con_LSF_1_LE_1./(100^2*1000),'ok','MarkerSize',2)
zlabel('Tau (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), LE: 1, (n = ',num2str(NumOfconLSF_1_LE_1_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100),...
    tau_min./(100^2*1000), tau_max./(100^2*1000)])

%con - LSF_1_LE_2
subplot(3,3,8)
plot3(alpha_con_LSF_1_LE_2,F_con_LSF_1_LE_2./(1000*100),...
    tau_con_LSF_1_LE_2./(100^2*1000),'ok','MarkerSize',2)
zlabel('Tau (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), LE: 2, (n = ',num2str(NumOfconLSF_1_LE_2_Files),')']);
axis([0, 2*pi, F_min./(1000*100), F_max./(1000*100),...
    tau_min./(100^2*1000), tau_max./(100^2*1000)])

set(fig7, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/LSF_1/Fig15_force_alpha_tau_LSF_1_LE_v12h.fig');
% saveas(gcf,'../Figures_MPC/LSF_1/Fig15_force_alpha_tau_LSF_1_LE_v12h.jpg');

disp('Figure 7 done')

%% Tau and phi on the same plot
% 
% figure(8);
% 
% load('CuratedData/LSF_1/mean_phi_con.mat')
% load('CuratedData/LSF_1/std_phi_con.mat')
% 
% %x-motion
% subplot(3,1,1)
% plot(Tstore,y_g,'-k','LineWidth',2)
% ylabel('y goal (cm)'); 
% xlabel('time (s)');
% axis([0, Tstore(end), -11, 11])
% 
% %con
% subplot(3,1,2)
% plot(time_hws,tau_con_LE_0./(100^2*1000),'ok','MarkerSize',2)
% hold on;
% plot(time_hws,zeroline,'-r','LineWidth',1)
% ylabel('Tau (N*m)'); 
% xlabel(['time (s) - control, (n = ',num2str(NumOfconLSF_p9_Files),')']);
% axis([0, Tstore(end), -0.01, 0.01])
% 
% %abdomen motion
% subplot(3,1,3)
% Fig1_19 = shadedErrorBar(Tstore,(180/pi)*mean_phi_con,(180/pi)*std_phi_con,'lineprops', '-k','transparent',false,'patchSaturation',0.075);
% Fig1_19.mainLine.LineWidth = 2;
% ylabel('phi (deg)');
% % plot(time_hws,mean_phi_con,'-k','Linewidth',2)
% hold on;
% plot(time_hws,zeroline,'-r','LineWidth',1)
% xlabel(['time (s) - control, (n = ',num2str(NumOfconLSF_p9_Files),')']);
% axis([0, Tstore(end), 200, 320])
% 
% saveas(gcf,'../Figures_MPC/LSF_1/Fig16_tau_phi_LSF_1_LE_v12h.fig');
% saveas(gcf,'../Figures_MPC/LSF_1/Fig16_tau_phi_LSF_1_LE_v12h.jpg');

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])