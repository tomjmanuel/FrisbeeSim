% iterate hyzer angle and visualize different lines (hopefully)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
format short
global m g Ia Id A d rho
global CLo CLa CDo CDa CMo CMa CRr
global CL_data CD_data CM_data CRr_rad CRr_AdvR CRr_data
global CMq CRp CNr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set "CoefUsed" = 1 OR 2
% This chooses the values of coefficients (specifies a set of if/then statements)
% to use for CLo CLa CDo CDa CMo CMa CRr.
% CoefUsed = 1 ... choose for using estimated short flights lift, drag, moment coefficients
% CoefUsed = 2 ... choose for using Potts and Crowther (2002) lift, drag, moment coefficients
CoefUsed=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define non-aerodynamic parameters
m = 0.175; % Kg
g = 9.7935; % m/s^2
A = 0.057; % m^2
d = 2*sqrt(A/pi); % diameter
rho = 1.23; % Kg/m^3
Ia = 0.002352; % moment of inertia about the spinning axis
Id = 0.001219; % moment of inertia about the planar axis'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THREE ESTIMATED COEFFICIENTS
%CMq= -0.005, CRp =-0.0055, CNr = 0.0000071 % short (three) flights
CMq= -1.44E-02; CRp =-1.25E-02; CNr = -3.41E-05; % long flight f2302
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE seven COEFFICIENTS estimated from three flights
CLo= 0.3331;
CLa= 1.9124;
CDo= 0.1769;
CDa= 0.685;
CMo= -0.0821;
CMa= 0.4338;
CRr= 0.00171; % for nondimensionalization = sqrt(d/g), magnitude of CRr changes
% with nondimensionalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THE seven COEFFICIENTS from Potts and Crowther (2002)
% CL =[ rad CL deg]

CL_data=[ -0.1745 -0.2250 -10;
-0.05236 0 -3;
0 0.150 0;
0.08727 0.4500 5;
0.17453 0.7250 10;
0.26180 0.9750 15;
0.34907 1.2000 20;
0.43633 1.4500 25;
0.52360 1.6750 30];

% CD =[ rad CD deg]
CD_data=[ -0.1745 0.1500 -10
-0.05236 0.0800 -3
0 0.1000 0
0.08727 0.1500 5
0.1745 0.2600 10
0.26180 0.3900 15
0.3491 0.5700 20
0.4363 0.7500 25
0.5236 0.9200 30];

% CM =[ rad CM deg]
CM_data=[-0.174532925 -0.0380 -10
-0.087266463 -0.0220 -5
-0.052359878 -0.0140 -3
0 -0.0060 0
0.052359878 -0.0060 3
0.104719755 -0.0020 6
0.157079633 0.0000 9
0.20943951 0.0100 12
0.261799388 0.0220 15
0.34906585 0.0440 20
0.401425728 0.0600 23
0.453785606 0.0840 26
0.523598776 0.1100 30];

% CRr_deg=[-5 -4 -3 -2 -1 0 1 2 3 4
% 5 6 7 8 9 10 11 12 13 14 15 30 ]
CRr_rad = [-0.0873 -0.0698 -0.0524 -0.0349 -0.0175 0.0000 0.0175 0.0349 0.0524 0.0698 0.0873 0.1047 0.1222 0.1396 0.1571 0.1745 0.1920 0.2094 0.2269 0.2443 0.2618 0.5236];
CRr_AdvR= [2 1.04 0.69 0.35 0.17 0];
CRr_data = [-0.0172 -0.0192 -0.018 -0.0192 -0.018 -0.0172 -0.0172 -0.0168 -0.0188 -0.0164 -0.0136 -0.01 -0.0104 -0.0108 -0.0084 -0.008 -0.008 -0.006 -0.0048 -0.0064 -0.008 -0.003 ...
-0.0112 -0.0132 -0.012 -0.0132 -0.012 -0.0112 -0.0112 -0.0108 -0.0128 -0.0104 -0.0096 -0.0068 -0.0072 -0.0076 -0.0052 -0.0048 -0.0048 -0.0028 -0.0032 -0.0048 -0.0064 -0.003 ...
-0.0056 -0.0064 -0.0064 -0.0068 -0.0064 -0.0064 -0.0052 -0.0064 -0.0028 -0.0028 -0.004 -0.002 -0.004 -0.002 -0.0016 0 0 0 0 -0.002 -0.0048 -0.003 ...
-0.0012 -0.0016 -0.0004 -0.0028 -0.0016 -0.0016 -0.0004 0.0004 0.0004 0.0008 0.0004 0.0008 0.0012 0.0008 0.002 0.0028 0.0032 0.0024 0.0028 0.0004 -0.0012 -0.003 ...
-0.0012 -0.0012 -0.0016 -0.0016 -0.0012 -0.0004 0.0004 0.0008 0.0008 0.0016 0.0004 0.002 0.0004 0.0016 0.002 0.002 0.002 0.0012 0.0012 0 -0.0012 -0.003 ...
-0.0012 -0.0012 -0.0004 -0.0008 -0.0008 -0.0008 0.0004 0.0004 0.0004 0.0008 0.0004 0.0008 -0.0004 0 0 0.0004 0 0 0.0004 -0.002 -0.0012 -0.003];

CRr_data = reshape(CRr_data,[22 6]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angle and angular velocities in rad and rad/sec respectively
% phi = angle about the x axis phidot = angular velocity
% theta = angle about the y axis thetadot = angular velocity
% gamma = angle about the z axis gd(gammadot) = angular velocity
% For reference, two sets of previously used initial conditions...
% Long flight (f2302) release conditions:
% thetao = 0.211; speedo = 13.5; betao = 0.15; gd=54
% Common release conditions:
% thetao = 0.192; speedo = 14; betao = 0.105; gd=50

% Define Simulation set initial conditions, enter your choosen values here:
thetao = .192; % initial pitch angle
speedo = 13.7; % magnitude, m/sec
betao = .105; % flight path angle in radians between velocity vector and horizontal
gd=50; % initial spin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphao = thetao - betao; % inital alpha
vxo = speedo * cos(betao); % velocity component in the nx direction
vzo = -(speedo * sin(betao)); % velocity component in the nz direction
% (note: nz is positive down)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0= vector of initial conditions
% x0= [ x y z vx vy vz phi theta phidot thetadot gd gamma]
% Specify the set of inital conditions to use:
% the first set of conditions is for a disc released:
% theta, speed, and spin as specified above (thetao, speedo, gd),
% 1 meter above the ground, forward and right 0.001,
% no roll angle, no velocity in the the y direction
% and for positive alpha, disc is pitched up, with a neg. w component
% the second set is the long flight f2302 estimated initial conditions 
% First set:
%x0= [ 0.001 0.001 -1 vxo 0 vzo 0 thetao 0.001 0.001 gd 0 ]
% Second set:
%x0= [ x y z vx vy vz phi theta phidot thetadot gd gamma]
%   x           y       z           vx         vy       vz      phi theta
%   

%%
roll_deg = -10:10; %degrees
nR = length(roll_deg);
nsteps = 300;
Paths = zeros(nR,2,nsteps+1); %[releaseangle x y]

for i=1:nR
    x0=[-9.03E-01 -6.33E-01 -9.13E-01 1.34E+01 -4.11E-01 1.12E-03 -7.11E-02 2.11E-01 -1.49E+01 -1.48E+00 5.43E+01 5.03E+00];
    x0(1:3)=0;

    %change "Roll angle"
    x0(7) = roll_deg(i)*(2*pi/360);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Enter values for tfinal and nsteps:
    tfinal = 1.46; % length of flight
    nsteps = 300; % number of time steps for data
    tspan=[0:tfinal/nsteps:tfinal];
    %options=[]
    %options = odeset('AbsTol', 0.000001,'RelTol', 0.00000001,'OutputFcn','odeplot');
    options = odeset('AbsTol', 0.000001,'RelTol', 0.00000001);
    % Calls the ODE and integrate the frisbee equations of motions in the
    % subroutine, discfltEOM.m
    [t,x]=ode45(@discfltEOM,tspan,x0,options,CoefUsed);
    xx = x(:,1);
    yy = x(:,2);
    Paths(i,1,:) = xx;
    Paths(i,2,:) = yy;
    
end

%% Visualize Paths
figure
for i=1:nR
    hold on
    plot(squeeze(Paths(i,1,:)),squeeze(Paths(i,2,:)))
    
end
