% The Folowing code is the demonstration for the labert transfer between
% 2 Planets of any choice with central body being Sun.
% Lambert code has been adopted based on Howard C. Curtis book algorithm
% pages 238-245 and pages 168-169
% Author: Adam Nekolny
%-------------------------------------------------------------------------
clear; clc; close all;
%-------------------------------------------------------------------------
% USER INPUT
origin_planet = 'Earth';
target_planet = 'Mars';
orbit_type = 'pro'; % pro - prograde, retro - retrograde
% Initial Time
year = 2036;
month = 7;
day = 6;
hours = 0;
minutes = 0;
% Arrival Date
year2 = 2037;
month2 = 5;
day2 = 24;
hours2 = 0;
minutes2 = 0;
% Animate Plot
anim = 0; % yes - 1 or no - 0
n = 10; % Number of spteps per animation frame
%-------------------------------------------------------------------------
% Segment Notification
fprintf('Running Lambert Transfer Calculation and Analysis\n');
fprintf('------------------------------------------------------------ \n');
% Constants
muS = 1.32712428e11; % km^3/s^2

% Load Ephemeris Origin Planet at Initial Time [km,km/s]
[r1i,v1i] = planetEphemeris(juliandate(year,month,day,hours,minutes,0),'Sun',origin_planet);
% Load Ephemeris Target Planet at Initial Time [km,km/s]
[r2i,v2i] = planetEphemeris(juliandate(year,month,day,hours,minutes,0),'Sun',target_planet);

% Load Ephemeris Origin Planet at Final Time [km,km/s]
[r1f,v1f] = planetEphemeris(juliandate(year2,month2,day2,hours2,minutes2,0),'Sun',origin_planet);
% Load Ephemeris Target Planet at Final Time [km,km/s]
[r2f,v2f] = planetEphemeris(juliandate(year2,month2,day2,hours2,minutes2,0),'Sun',target_planet);

% Time of Flight [s]
tof = 24*3600*(juliandate(year2,month2,day2,hours2,minutes2,0) - juliandate(year,month,day,hours,minutes,0));

% Lamber Transfer
[V1, V2] = lambert(r1i,r2f,tof,orbit_type,muS); % 

% Setting Up Planets as Initial Conditions
ICP1 = [r1i,v1i]; %[x,y,z,vx,vy,vz]
ICP2 = [r2i,v2i];
ICS1 = [r1i,V1];

% Setting up ODE45
options = odeset('RelTol',1e-13,'AbsTol',1e-13);

% Integration Time
ti = 0;
dt = 3600; % 1 hour
time = ti:dt:tof;

% Propagation fo Planets and Spacecraft In Parallel
% Note: This is to speed up the integration through using simultaneous
% computation in paralel, which works only for independent data sets. 
f1 = parfeval(@ode45, 2, @dR_2body, time, ICP1,options);
f2 = parfeval(@ode45, 2, @dR_2body, time, ICP2,options);
f3 = parfeval(@ode45, 2, @dR_2body, time, ICS1,options);

[t1, RVP1] = fetchOutputs(f1); % Results from the first instance
[t2, RVP2] = fetchOutputs(f2); % Results from the second instance
[t3, RVS1] = fetchOutputs(f3); % Results from the third instance

% Here I would have integration of the planets between the Nov 7 and Nov 13
% tof_wait = 
% f1 = parfeval(@ode45, 2, @dR_2body, time_wait, RVP1(end,:),options);
% f2 = parfeval(@ode45, 2, @dR_2body, time_wait, RVP2(end,:),options);

% 
% % Determine the ner time of flight
% tof2 = 24*3600*(juliandate(2036,7,6,0,0,0) - juliandate(year2,month2,day2,hours2,minutes2,0));
% time2 = ti:dt:tof2;
% 
% % Second Integration just planets to get the new planet positions
% f4 = parfeval(@ode45, 2, @dR_2body, time2, RVP1(end,:),options);
% f5 = parfeval(@ode45, 2, @dR_2body, time2, RVP2(end,:),options);
% 
% [t1, RVP12] = fetchOutputs(f4); % Results from the first instance
% [t2, RVP22] = fetchOutputs(f5); % Results from the second instance
% 
% 
% % Then perform labert or hohmann for the spacecraft for interplnaetary
% % transfer
% % This would be either new lambert or new hohman from mars back to earth
% [V12, V22] = lambert(RVP2(end,1:3),RVP12(end,1:3),tof2,orbit_type,muS); % 
% 
% ICS12 = [RVP2(end,1:3),V12];
% % Integrate for the spacecrfat
% f6 = parfeval(@ode45, 2, @dR_2body, time2, ICS12,options);
% [t3, RVS12] = fetchOutputs(f6); % Results from the third instance

if anim == 1
    % Initialize figure
    figure;
    hold on;
    axis equal;
    view(3);
    
    % Create Sun as a sphere at the origin
    [sx, sy, sz] = sphere(50); % Create sphere coordinates
    sun_radius = 6.96e6*2; % Approximate radius of the Sun in meters
    sun = surf(sx * sun_radius, sy * sun_radius, sz * sun_radius, ...
        'EdgeColor', 'none', 'FaceColor', 'yellow', 'FaceLighting', 'gouraud', 'DisplayName', 'Sun');
    
    % Initialize plot lines for the origin planet, target planet, and Lambert transfer
    origin_line = plot3(NaN, NaN, NaN, 'LineWidth', 2, 'DisplayName', origin_planet);
    target_line = plot3(NaN, NaN, NaN, 'LineWidth', 2, 'DisplayName', target_planet);
    transfer_line = plot3(NaN, NaN, NaN, 'LineWidth', 2, 'DisplayName', 'Lambert Transfer');
    
    % Set labels and title
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Orbit Animation With Lambert Transfer');
    grid on;
    legend show;
    
    % Determine the number of time steps
    num_steps = size(RVP1, 1);
    % Loop over time steps to animate the trajectories
    for k = 1:n:num_steps
        % Update the data for each line up to the current time step
        set(origin_line, 'XData', RVP1(1:k, 1), 'YData', RVP1(1:k, 2), 'ZData', RVP1(1:k, 3));
        set(target_line, 'XData', RVP2(1:k, 1), 'YData', RVP2(1:k, 2), 'ZData', RVP2(1:k, 3));
        set(transfer_line, 'XData', RVS1(1:k, 1), 'YData', RVS1(1:k, 2), 'ZData', RVS1(1:k, 3));
        
        % Refresh the plot
        drawnow;
        
        % Optional: Pause to control the speed of the animation
    end
    
    hold off;
else
    % Plot Trajectories and add the Sun as a sphere at the origin 3D
    figure;
    hold on
    axis equal;
    view(3);
    [sx, sy, sz] = sphere(50); % Create sphere coordinates
    sun_radius = 6.96e6; % Approximate radius of the Sun in meters (adjust scale if needed)
    surf(sx * sun_radius, sy * sun_radius, sz * sun_radius,'EdgeColor', 'none', 'FaceColor', 'yellow', 'FaceLighting', 'gouraud', 'DisplayName','Sun'); % Sun visualization
    plot3(RVP1(:,1),RVP1(:,2),RVP1(:,3), 'LineWidth', 2, 'DisplayName',origin_planet);
    %plot3(RVP12(:,1),RVP12(:,2),RVP12(:,3), 'LineWidth', 2, 'DisplayName',origin_planet);
    plot3(RVP2(:,1),RVP2(:,2),RVP2(:,3), 'LineWidth', 2, 'DisplayName',target_planet);
   % plot3(RVP22(:,1),RVP22(:,2),RVP22(:,3), 'LineWidth', 2, 'DisplayName',target_planet);
    plot3(RVS1(:,1),RVS1(:,2),RVS1(:,3), 'LineWidth', 2, 'DisplayName','Lambert Transfer');
    %plot3(RVS12(:,1),RVS12(:,2),RVS12(:,3), 'LineWidth', 2, 'DisplayName','Lambert Transfer');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Orbit Plot With Lambert Transfer');
    grid on
    legend show;
    hold off
end

% Plot Trajectories in 2D
figure;
hold on;
axis equal
tht = linspace(0, 2*pi, 50);
fill(6.96e6*cos(linspace(0, 2*pi, 50)), 6.96e6*sin(linspace(0, 2*pi, 50)), 'y','EdgeColor','y','DisplayName','Sun');
plot(RVP1(:,1),RVP1(:,2),'LineWidth',2, 'DisplayName',origin_planet);
plot(RVP2(:,1),RVP2(:,2),'LineWidth',2, 'DisplayName',target_planet);
plot(RVS1(:,1),RVS1(:,2),'LineWidth',2, 'DisplayName','Lambert Transfer');
xlabel('X'); ylabel('Y');
title('2D Orbit Plot With Lambert Tranfer');
legend show;
hold off;

% ------------------------------------------------------------------------
% Delta-V Calculation of Lamber transfer
fprintf('Delta-V Lambert\n');
fprintf('Note: Assuming Launch from SOI \n');
% Lamber Delta-V Calculation
nVP1 = norm(v1i); 
nVS1 = norm(V1);
delV1 = nVS1-nVP1; % Delta-V required at original planet

nVS2 = norm(V2);
RVP1E = RVP1(end,4:6);
nVP2 = norm(RVP1E);
delV2 = nVP2-nVS2; % delta-v required at target planet

tot_delVLam = abs(delV2)+abs(delV1); % total delta-v
fprintf('\n');
fprintf('Delta-V Original Planet : %.4f km/s\n',delV1);
fprintf('Delta-V Target Planet : %.4f km/s\n',delV2);
fprintf('Total Delta-V : %.4f km/s\n',tot_delVLam);
fprintf('\n');

% Time of flight
fprintf('Time of Flight Lambert\n');
fprintf('TOF: %0.2f Days\n',tof/(24*3600));
fprintf('------------------------------------------------------------ \n');
% ------------------------------------------------------------------------
% Error Analysis
fprintf('Comparing Position and Velocity Integration to Ephemeris \n');

% Relative errors for final position and velocity of planets
rp1_rel_err = abs(norm(r1f)-norm(RVP1(end,1:3)))/norm(r1f) * 100;
rp2_rel_err = abs(norm(r2f)-norm(RVP2(end,1:3)))/norm(r2f) * 100;
vp1_rel_err = abs(norm(v1f)-norm(RVP1(end,4:6)))/norm(v1f) * 100;
vp2_rel_err = abs(norm(v2f)-norm(RVP2(end,4:6)))/norm(v2f) * 100;

fprintf('%s \n',origin_planet);
fprintf('Error with Respect to its Actual Ephemeris is: \n');
fprintf('Relative Position Error: %.4f %% \n',rp1_rel_err);
fprintf('Relative Velocity Error: %.4f %% \n',vp1_rel_err);
fprintf('\n');

fprintf('%s \n',target_planet);
fprintf('Error for Target Planet with Respect to its Actual Ephemeris is: \n');
fprintf('Relative Position Error: %.4f %% \n',rp2_rel_err);
fprintf('Relative Velocity Error: %.4f %% \n',vp2_rel_err);
fprintf('\n');


fprintf('------------------------------------------------------------ \n');
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% SUB-FUNCTIONS USED WITHIN THE CODE
%-----------------------------------------------------------------------
% Howard Curtis Orbital Mechanics For Engineering Students, D25
function [V1, V2] = lambert(R1, R2, t, string,mu) 
%-----------------------------------------------------------------------
%{
    This function solves Lambert’s problem.
    mu - gravitational parameter (km^3/s^2)
    R1, R2 - initial and final position vectors (km)
    r1, r2 - magnitudes of R1 and R2
    t - the time of flight from R1 to R2 (a constant) (s)
    V1, V2 - initial and final velocity vectors (km/s)
    c12 - cross product of R1 into R2
    theta - angle between R1 and R2
    string - ’pro’ if the orbit is prograde
    ’retro’ if the orbit is retrograde
    A - a constant given by Equation 5.35
    z - alpha*x^2, where alpha is the reciprocal of the
    semimajor axis and x is the universal anomaly
    y(z) - a function of z given by Equation 5.38
    F(z,t) - a function of the variable z and constant t,
    - given by Equation 5.40
    dFdz(z) - the derivative of F(z,t), given by Equation 5.43
    ratio - F/dFdz
    tol - tolerance on precision of convergence
    nmax - maximum number of iterations of Newton’s procedure
    f, g - Lagrange coefficients
    gdot - time derivative of g
    C(z), S(z) - Stumpff functions
    dum - a dummy variable
    User M-functions required: stumpC and stumpS
%}
% ––––––––––––––––––––––––––––––––––––––––––––––
%...Magnitudes of R1 and R2:
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1, R2);
theta = acos(dot(R1,R2)/r1/r2);

%...Determine whether the orbit is prograde or retrograde:
if nargin < 4 jj (~strcmp(string,'retro') & (~strcmp(string,'pro')))
    string = 'pro';
    fprintf('\n ** Prograde trajectory assumed.\n')
end

if strcmp(string,'pro')
    if c12(3) <= 0
        theta = 2*pi - theta;
    end
elseif strcmp(string,'retro')
    if c12(3) >= 0
        theta = 2*pi - theta;
    end
end

%...Equation 5.35:
A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));

%...Determine approximately where F(z,t) changes sign, and
%...use that value of z as the starting value for Equation 5.45:
z = -100;
while F(z,t) < 0
    z = z + 0.1;
end

%...Set an error tolerance and a limit on the number of iterations:
tol = 1.e-8;
nmax = 5000;

%...Iterate on Equation 5.45 until z is determined to within the
%...error tolerance:
ratio = 1;
n = 0;
while (abs(ratio) > tol) & (n <= nmax)
    n = n + 1;
    ratio = F(z,t)/dFdz(z);
    z = z - ratio;
end

%...Report if the maximum number of iterations is exceeded:
if n >= nmax
    fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax)
end

%...Equation 5.46a:
f = 1 - y(z)/r1;

%...Equation 5.46b:
g = A*sqrt(y(z)/mu);

%...Equation 5.46d:
gdot = 1 - y(z)/r2;

%...Equation 5.28:
V1 = 1/g*(R2 - f*R1);

%...Equation 5.29:
V2 = 1/g*(gdot*R2 - R1);

return
%-----------------------------------------------------------------------
% Subfunctions used in the main body:
%-----------------------------------------------------------------------
%...Equation 5.38:
function dum = y(z)
    dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
end

%...Equation 5.40:
function dum = F(z,t)
    dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
end

%...Equation 5.43:
function dum = dFdz(z)
    if z == 0
        dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
    else
        dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) + A*sqrt(C(z)/y(z)));
    end
end

%...Stumpff functions Equations 3.52-3.53:
function dum = C(z)
    if z > 0
        dum = (1-cos(sqrt(z)))/z;
    elseif z < 0
        dum = (cosh(sqrt(-z))-1)/(-z);
    else
        dum = 1/2;
    end
end

function dum = S(z)
    if z > 0
        dum = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        dum = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z))^3;
    else
        dum = 1/6;
    end
end

end %lambert
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% 2 Body Integrator from twobody_integrator by Erik Butcher Ph.D.
function dRV = dR_2body(time,RV)

x  = RV(1);
y  = RV(2);
z  = RV(3);
vx = RV(4);
vy = RV(5);
vz = RV(6);
r=sqrt(x^2+y^2+z^2);

mu=1.32712428e11; % [km^3/s^2]

ax = -mu*x/r^3;
ay = -mu*y/r^3;
az = -mu*z/r^3;

dRV = [vx; vy; vz; ax; ay; az];

end
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% END
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx






