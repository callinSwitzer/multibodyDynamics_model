%4/4/17- Setting up my model (Model 1) in MATLAB using ode45
    %Models 2-7 since then (and even more errors/fixes)
    %Model 8 was essentially refined on 6/29/17
%7/6/17- Model 9 development will incorporate our cost function and 
%multiple wing strokes.
%8/7/17- Model 10(a?) will do all of the above for SEVERAL full run 
%throughs. Model 10(b?) will have a Gaussian cost function. a.k.a. no more
%window!
    %9/4/17- FIXED K and c
    %9/6/17- Fixed tau to have a maximum value of 10000 g*(cm^2)/(s^2)
    %2/8/18- Eliminated the need for global variables and implemented 
        %structs that are passed into the ODE.
%3/6/18- Model 11 is born.
        %Restructured the ENTIRE layout of the variables to have nested
        %structs which contain what we want.
        %hws contains a 2500x1 struct of: 
            %ICs, F, alpha, tau, Q, cost, NewICs
        %Each hws is nested in a 100x1 large struct for Qstore and Winstore
    %7/3/18- Fixed a MAJOR error in the resting angle of the torsional spring.
        %Now there is the appropriate restorative force.
%7/18/18- Model 12 was born.
        %Fixed another MAJOR error. tau0 is now a fixed force (as
        %originally intended), needed to fix myODE and remove the
        %sinusoidal component of the applied torque.
    %~8/17/18- Incorporated LengthScaleFactor to change the size of the model.
    %10/16/18- Incorporated Size
    %1/7/19- Finally attempting a formal Model Predictive Control (MPC)
    %1/31/19- Corrected issues pertaining to the Moment of Inertia.
        %myODE5 was born.
        %1) incorrect multiplication in myODE4,
        %2) accounting for parallel axis theorem (due to petiole length
        %extension) in myODE4.
    %2/14/19- For the sake of efficiency, lines 184, and 314 were muted
        %to no longer create the Qstore struct. 

%12b is for horizontal aggressive maneuver
%12c is for vertical aggressive maneuver
%12h is for sum of prime number sines

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
start_time = datetime;
disp(['Time at start: ',datestr(start_time)])
tic

%% Variable setup

LengthScaleFactor = 1; %This will multiply all linear scales of the model 
    %appropriately to modify the size of the model
LengthExtend = 0; %This will lengthen the distance between the two masses

bl = 2*LengthScaleFactor*(0.908+1.7475); %This is the bodylength of the 
           %moth in cm
           
%Pack the parameter struct (NOT to be altered)
par.L1 = LengthScaleFactor*0.908; %Length from the thorax-petiole joint to 
    %the center of the head-thorax in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.47195), BEE (0.2915), MOTH SCENARIO (0.908)
par.L_petiole = LengthExtend*bl; %Length of petiole extension as a 
    %percentage of body length.
par.L2 = LengthScaleFactor*1.7475 + par.L_petiole; 
    %Length from the thorax-petiole joint 
    %to the center of the gaster in cm - THIS WILL CHANGE 
          %FOR DAUBER (1.177466667), BEE (0.2935), MOTH SCENARIO (1.7475)
par.L3 = LengthScaleFactor*0.75; %Length from the thorax-petiole joint to 
    %the aerodynamic force vector in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.25), BEE (0.25), MOTH SCENARIO (0.75)
par.rho = 1; %The density of the insect in g/(cm^3)
par.rhoA = 1.18*10^-3; %The density of air in g/(cm^3)
par.muA = 1.86*10^-4; %The dynamic viscosity of air at 27C in g/(cm*s)
par.ahead = LengthScaleFactor*0.908; %Major axis of head-thorax ellipsoid 
          %in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.47195), BEE (0.2915), MOTH SCENARIO (0.908)
par.abutt = LengthScaleFactor*1.7475; %Major axis of abdomen ellipsoid
          %in cm  - THIS WILL CHANGE 
          %FOR DAUBER (0.35), BEE (0.2935), MOTH SCENARIO (1.7475)
par.bhead = LengthScaleFactor*0.507; %Minor axis of head-thorax ellipsoid 
          %in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.1839665), BEE (0.1815), MOTH SCENARIO (0.507)
par.bbutt = LengthScaleFactor*0.1295; %Minor axis of abdomen ellipsoid 
          %in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.15), BEE (0.173), MOTH SCENARIO (0.1295)
par.K = LengthScaleFactor*29.3; %K is the torsional spring constant of 
          %the thorax-petiole joint in (cm^2)*g/(rad*(s^2))
par.c = LengthScaleFactor*14075.8; %c is the torsional damping constant of 
            %the thorax-petiole joint in (cm^2)*g/s
par.g = 980; %g is the acceleration due to gravity in cm/(s^2)

%% More variables NOT to alter
halfwingStrokes = 501; %Originally 100, but we need an extra half wing 
          %stroke for the formal MPC. 
hwbf = 50; %wing beat frequency of the insect in Hz - THIS WILL CHANGE 
          %FOR DAUBER (100), BEE (400), MOTH SCENARIO (50)
timestep = 100; %The number of timesteps we will use when interpolating
                %DO NOT CHANGE TIMESTEP
kk = 0; %Overall counter. Do not reset!
i_overall = 0; %Overall counter of runs. DO NOT RESET!
kkskip = 1; %Overall counter OF SKIPS. Do not reset!
nn = 1; %Counter of full runs. Do not reset!
MaxWorkers = feature('numcores'); 
    %Number of workers for parallel computing. Do not reset!
RecedeFrac = 0.25; %Must be a value between 0 and 1.
                  %This governs how far back the goal is shifted from the 
                  %original goal.   
goal_index = timestep:(RecedeFrac*100):(timestep*halfwingStrokes); 
    %Indices for the goal (without having to calculate it every loop)
    
%% Variables to alter (if someone else is making changes for me)
numOfTrajectories = 2500; %Number of trajectories per half wingstroke spray
shift = 0; %This is to increase the file number as appropriate
FullRuns = 10; %Number of full runs 
                  
%% Time vector
Tstore = linspace(0,((1/hwbf)*halfwingStrokes),...
    (timestep*halfwingStrokes))';
% save('12h-AMSumOfPrimes_MPC/Tstore_MPC_hws_sp.mat','Tstore');

%% Relevant to cost function
%Goal criteria
x_g = 0; %in cm
theta_g = pi/4; %in rad

xdot_g = 0; %in cm/s
thetadot_g = 0; %in rad/s

%Create the sum of primes signal
signal_amp = 5; %in cm THIS WAS 5! I ALTERED THIS WHEN SETTING UP ARDUINO
prime_f = [0.2, 0.3, 0.5, 0.7, 1.1, 1.7, 2.9, 4.3, 7.9, 13.7, 19.9]; %in Hz
prime_a = (signal_amp./(2*pi.*prime_f)).*(2*pi.*prime_f(1)); %in cm
prime_ph = zeros(1,numel(prime_f)); %Not sure if we'll require a phase

%y-motion goal criteria (for the cost function)
y_g = prime_a(1)*sin(2*pi*prime_f(1)*Tstore + prime_ph(1)) +...
    prime_a(2)*sin(2*pi*prime_f(2)*Tstore + prime_ph(2)) +...
    prime_a(3)*sin(2*pi*prime_f(3)*Tstore + prime_ph(3)) +...
    prime_a(4)*sin(2*pi*prime_f(4)*Tstore + prime_ph(4)) +...
    prime_a(5)*sin(2*pi*prime_f(5)*Tstore + prime_ph(5)) +...
    prime_a(6)*sin(2*pi*prime_f(6)*Tstore + prime_ph(6)) +...
    prime_a(7)*sin(2*pi*prime_f(7)*Tstore + prime_ph(7)) +...
    prime_a(8)*sin(2*pi*prime_f(8)*Tstore + prime_ph(8)) +...
    prime_a(9)*sin(2*pi*prime_f(9)*Tstore + prime_ph(9)) +...
    prime_a(10)*sin(2*pi*prime_f(10)*Tstore + prime_ph(10)) +...
    prime_a(11)*sin(2*pi*prime_f(11)*Tstore + prime_ph(11)); %in cm

ydot_g = 2*pi*prime_a(1)*prime_f(1)*cos(2*pi*prime_f(1)*Tstore + prime_ph(1)) +...
    2*pi*prime_a(2)*prime_f(2)*cos(2*pi*prime_f(2)*Tstore + prime_ph(2)) +...
    2*pi*prime_a(3)*prime_f(3)*cos(2*pi*prime_f(3)*Tstore + prime_ph(3)) +...
    2*pi*prime_a(4)*prime_f(4)*cos(2*pi*prime_f(4)*Tstore + prime_ph(4)) +...
    2*pi*prime_a(5)*prime_f(5)*cos(2*pi*prime_f(5)*Tstore + prime_ph(5)) +...
    2*pi*prime_a(6)*prime_f(6)*cos(2*pi*prime_f(6)*Tstore + prime_ph(6)) +...
    2*pi*prime_a(7)*prime_f(7)*cos(2*pi*prime_f(7)*Tstore + prime_ph(7)) +...
    2*pi*prime_a(8)*prime_f(8)*cos(2*pi*prime_f(8)*Tstore + prime_ph(8)) +...
    2*pi*prime_a(9)*prime_f(9)*cos(2*pi*prime_f(9)*Tstore + prime_ph(9)) +...
    2*pi*prime_a(10)*prime_f(10)*cos(2*pi*prime_f(10)*Tstore + prime_ph(10)) +...
    2*pi*prime_a(11)*prime_f(11)*cos(2*pi*prime_f(11)*Tstore + prime_ph(11)); %in cm/s

%Weighting coefficients
%c1 = xdot, c2 = ydot, c3 = thetadot, c4 = x, c5 = y, c6 = theta
% c1 = 1*10^-5; c2 = 1*10^-5; c3 = 10^6; c4 = 10^7; c5 = 10^8; c6 = 10^10; 
%New weighting functions (as of 2/3/19)
c1 = 1*10^-5; c2 = 1*10^-5; c3 = 10^8; c4 = 10^9; c5 = 10^10; c6 = 10^10; 

%% Initial conditions
  %q0 = [x, y, theta,   phi,            xdot,    ydot,    thetadot, phidot]
q0_og = [0; 0; theta_g; (theta_g + pi); 1*10^-4; 1*10^-4; 0; 0];
q0 = q0_og;
%Options for ode45
OPTIONS = odeset('RelTol',1e-2,'AbsTol',1e-4); 

par.betaR = q0_og(4) - q0_og(3) - pi; %This is the resting configuration of our 
    %torsional spring(s) = Initial abdomen angle - initial head angle - pi

%% Prescribe the storage data
% xx = zeros(timestep*halfwingStrokes,8);
xx = zeros(timestep,numel(q0_og));
% Q = [];
tt = linspace(0,1/hwbf,timestep)';
for i = 1:((halfwingStrokes-1)*(1/RecedeFrac))
    strucnames{i,:}= ['PartPath',num2str(i)];
end

for i = 1:((halfwingStrokes-1)*(1/RecedeFrac))
%     Qstore{i,1} = struct(strucnames{i},[]);
    Winstore{i,1} = struct(strucnames{i},[]);
end

cost_perspray = zeros(numOfTrajectories,1);

timing1 = toc;

%% Using ode45 
tic
%Start the full run loop
for nn = 1:FullRuns

%Start the iteration of consecutive wing strokes
for i = 1:((halfwingStrokes-1)*(1/RecedeFrac))
    %Overall counter
    i_overall = i_overall+1;
    
%To estimate how long the run will take
%First timestamp
    if i_overall == 1
        timeMarker(1) = datetime;
    end
%First estimate: 1 partial path
    if i_overall == 2
        timeMarker(2) = datetime;
        timeDiff(1) = timeMarker(2)-timeMarker(1);
        timeEst(1) = timeDiff(1)*((halfwingStrokes-1)*...
            (1/RecedeFrac))*FullRuns + start_time;
    end 
%Second estimate: 100 partial paths
    if i_overall == 100
        timeMarker(3) = datetime;
        timeDiff(2) = timeMarker(3)-timeMarker(1);
        timeEst(2) = timeDiff(2)*((halfwingStrokes-1)*...
            (1/RecedeFrac)/100)*FullRuns + start_time;
    end 
%Third estimate: 1/10th of the way for all full runs
    if i_overall == ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns/10)
        timeMarker(4) = datetime;
        timeDiff(3) = timeMarker(4)-timeMarker(1);
        timeEst(3) = timeDiff(3)*10 + start_time;
    end
%Fourth estimate: Halfway there (for all full runs)
    if i_overall == ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns/2)
        timeMarker(5) = datetime;
        timeDiff(4) = timeMarker(5)-timeMarker(1);
        timeEst(4) = timeDiff(4)*2 + start_time;
    end
%Fifth (and final) estimate: 3/4 of the way there (for all full runs)
    if i_overall == ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns*0.75)
        timeMarker(6) = datetime;
        timeDiff(5) = timeMarker(6)-timeMarker(1);
        timeEst(5) = timeDiff(5) + timeDiff(5)/3 + start_time;
            %The logic here is: 75% + (75/3)%
    end
    
rng('shuffle'); %This re-seeds the random number generator

PartPath(1:numOfTrajectories) = struct('ICs',q0,'F',0,'alpha',0,'tau0',0,...
    'bigQ',[],'cost',NaN,'NewICs',[],'check',[]);

    for r_i = 1:numOfTrajectories
        PartPath(r_i).F = 44300*rand(1)*(LengthScaleFactor^3); %For LSF >= 1
%         PartPath(r_i).F = 44300*rand(1); %For LSF < 1
            %F is the aerodynamic vector in g*cm/(s^2) 44300 for hawkmoth
        PartPath(r_i).alpha = 2*pi*rand(1); %alpha is the angle of the 
            %aerodynamic vector with respect to the head-thorax mid-line 
            %in radians
        PartPath(r_i).tau0 = 100000*10*(2*(rand(1)-0.5))*(LengthScaleFactor^4); %For LSF >= 1
%         PartPath(r_i).tau0 = 100000*(2*(rand(1)-0.5)); %For LSF < 1
            %The initial torque 
            %applied in g*(cm^2)/(s^2) old var -> 10*10000
%         PartPath(r_i).tau0 = 0; %The initial torque applied in g*(cm^2)/(s^2)
    end
    
%Define the goal criteria for y and ydot BEFORE the parfor loop
ydot_g_thisPartPath = ydot_g(goal_index(i),1);
y_g_thisPartPath = y_g(goal_index(i),1);
    
    parfor(j = 1:numOfTrajectories, MaxWorkers)   
       %[TOUT,YOUT,TE,YE,IE] = ode45(ODEFUN,TSPAN,Y0,OPTIONS)
        [T, Q] = ode45(@(t,Q) myODE_5(t, Q, par,...
            PartPath(j).F, PartPath(j).alpha, PartPath(j).tau0),...
            [0, 1/hwbf], PartPath(j).ICs, OPTIONS);
        
        xx = interp1(T, Q, tt);
        bigQ{j,1} = xx;

        %Cost function
        cost(j,1) = c1*(xx(end,5) - xdot_g)^2 +...
            c2*(xx(end,6) - ydot_g_thisPartPath)^2 +...
            c3*(xx(end,7) - thetadot_g)^2 +...
            c4*(xx(end,1) - x_g)^2 +...
            c5*(xx(end,2) - y_g_thisPartPath)^2 +...
            c6*(xx(end,3) - theta_g)^2;
    
        %Store the control, state vars, and cost which satisfy the conditions
        %Keep in mind, the new initial conditions will now be shifted
        %backward to the receding horizon fraction. (Cost function will
        %still determine which path to take).
        NewICs{j,1} = [xx((end*RecedeFrac+1),1); xx((end*RecedeFrac+1),2);...
            xx((end*RecedeFrac+1),3); xx((end*RecedeFrac+1),4); ...
            xx((end*RecedeFrac+1),5); xx((end*RecedeFrac+1),6); ...
            xx((end*RecedeFrac+1),7); xx((end*RecedeFrac+1),8)];
        %NOTE: Index is offset an additional one so that each signal
            %extracted later is 25 counts (as opposed to 24). 
            %Because of the repeated values from the last state
            %variables and the first of the next partial path, we must
            %ignore the first value of each partial path when we
            %extract the values to create a continuous signal later.
        
        check{j,1} = [PartPath(j).F, PartPath(j).alpha, PartPath(j).tau0];

    end %end the parfor

kk = kk+numOfTrajectories; %Overall counter (in bulk). DO NOT RESET
    
%To clear the pesky warnings/errors
clc;
        
    for a_i = 1:numOfTrajectories
        PartPath(a_i).bigQ = bigQ{a_i,1};  %This is the set of all trajectories
        PartPath(a_i).cost = cost(a_i,1);  %A column devoted to cost
        PartPath(a_i).NewICs = NewICs{a_i,1};  %The new initial conditions
        PartPath(a_i).check = check{a_i,1}; %This is to check that the applied 
            %efforts are in the proper order.
    end
    
%Store the set of realizations into Qstore
% Qstore{i,1} = PartPath;
    
    %Unpack cost values to identify the lowest cost. A separate loop is
    %necessary to do this step.
    for c_i = 1:numOfTrajectories
        cost_perspray(c_i,1) = PartPath(c_i).cost;
    end

%Display the progress of the program
disp(['Time at start: ',datestr(start_time)])
disp('Simulation in progress...')
disp(['Full run ', num2str(nn), ' of ', num2str(FullRuns)]);
disp(['Each half wingstroke is divided into ', num2str(1/RecedeFrac), ' partial paths.']);
disp(['Partial path ', num2str(i), ' of ', num2str((halfwingStrokes-1)*(1/RecedeFrac))]);
disp(['Half wingstroke ', num2str(i*RecedeFrac), ' of ', num2str(halfwingStrokes-1)]);
disp(' ');
if i_overall >= 2 && i_overall < 100
    disp(['Estimated end time of simulation is: ', datestr(timeEst(1))]);
end

if i_overall >= 100 && i_overall < ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns/10)
    disp(['Estimated end time of simulation is: ', datestr(timeEst(2))]);
end

if i_overall >= ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns/10) &&...
        i_overall < ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns/2)
    disp(['Estimated end time of simulation is: ', datestr(timeEst(3))]);
end

if i_overall >= ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns/2) &&...
        i_overall < ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns*0.75)
    disp(['Estimated end time of simulation is: ', datestr(timeEst(4))]);
end

if i_overall >= ((halfwingStrokes-1)*(1/RecedeFrac)*FullRuns*0.75)
    disp(['Estimated end time of simulation is: ', datestr(timeEst(5))]);
end

disp(' ');
disp([num2str(kkskip-1), ' loops have been skipped']);
disp(['You are ',... 
    num2str(100*kk/(numOfTrajectories*(halfwingStrokes-1)*(1/RecedeFrac)*FullRuns)),...
    '% of the way there!']);

    if kk == numOfTrajectories*halfwingStrokes*FullRuns
        disp('All trajectories are complete. No instabilities!')
    end

LowestCost_index = find(cost_perspray==min(cost_perspray));

    if numel(LowestCost_index) < 1
        kkskip = kkskip + 1;
        break 
    end
    
%Cost index identifies the winner of this set of realizations. Thus, the
%values stored in Winstore are as follows

Winstore{i,1} = PartPath(LowestCost_index);

%Our NEW initial conditions in the following order:
%Q = [x, y, theta,   phi,  xdot,    ydot, thetadot, phidot]
q0(1:8,1) = PartPath(LowestCost_index).NewICs;

%In the event that q0 ends up being all NaNs
    if isnan(q0) == 1
        kkskip = kkskip + 1;
        break 
    end
    
end %End of loop for the current half wing stroke

%SAVE the vars necessary for big data stuff
%Six suffixes relevant for our tests:
%sp = aggressive maneauver, Sum of Prime frequencies
%con = control
%ntc = control, no tau
%nts = no tau, stiff
%ntf = no tau, floppy
%yts = yes tau, stiff
%ytf = yes tau, floppy

% save(['12h-AMSumOfPrimes_MPC/Qstore_',num2str(nn+shift),'_hws_sp_con.mat'],'Qstore');
if LengthScaleFactor >= 1
    if LengthExtend < 1 && LengthExtend > 0
    save(['12h-AMSumOfPrimes_MPC/Winstore_',num2str(nn+shift),...
        '_hws_sp_con_LSF_',num2str(LengthScaleFactor),...
        '_LE_p',num2str(10*LengthExtend),'.mat'],'Winstore');
    else
    save(['12h-AMSumOfPrimes_MPC/Winstore_',num2str(nn+shift),...
        '_hws_sp_con_LSF_',num2str(LengthScaleFactor),...
        '_LE_',num2str(LengthExtend),'.mat'],'Winstore');
    end
else
    if LengthExtend < 1 && LengthExtend > 0
    save(['12h-AMSumOfPrimes_MPC/Winstore_',num2str(nn+shift),...
        '_hws_sp_con_LSF_p',num2str(10*LengthScaleFactor),...
        '_LE_p',num2str(10*LengthExtend),'.mat'],'Winstore');
    else
    save(['12h-AMSumOfPrimes_MPC/Winstore_',num2str(nn+shift),...
        '_hws_sp_con_LSF_p',num2str(10*LengthScaleFactor),...
        '_LE_',num2str(LengthExtend),'.mat'],'Winstore');
    end
end
%Now we reset for the next full run
%Our initial conditions in the following order:
%Q = [x, y, theta,   phi,  xdot,    ydot, thetadot, phidot]
q0 = q0_og;

%Prescribe the storage data
xx = zeros(timestep*halfwingStrokes,8);

end %Finished will ALL runs.

timing2 = toc;

%% Figures of state variables, trajectories, and flexion w.r.t. time.
disp(['Loading vars took ', num2str(timing1), ' seconds'])

if timing2>60 && timing2<180
    disp([num2str(FullRuns),' run(s) took ', num2str(timing2), ' seconds'])
elseif timing2>=180 && timing2<3600
    disp([num2str(FullRuns),' run(s) took ', num2str(timing2/60), ' minutes'])
elseif timing2>=3600 && timing2<72*3600
    disp([num2str(FullRuns),' run(s) took ', num2str(timing2/(60*60)), ' hours'])
elseif timing2>=72*3600
    disp([num2str(FullRuns),' run(s) took ', num2str(timing2/(24*60*60)), ' days'])
end

%% Final display on window
disp(' ');
disp('Estimated end time of simulations were as follows... ')
disp(['First estimate (One partial path):         ', datestr(timeEst(1))]);
disp(['Second estimate (First 100 partial paths): ', datestr(timeEst(2))]);
disp(['Third estimate (1/10th of the way done):   ', datestr(timeEst(3))]);
disp(['Fourth estimate (Halfway done):            ', datestr(timeEst(4))]);
disp(['Fifth estimate (3/4 of the way done):      ', datestr(timeEst(5))]);
disp(' ');
disp('Completely done, son.');
disp(['Time at end: ',datestr(datetime)])