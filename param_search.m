% Adjust disc stability and visualize lines hopefully
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
format short
global m g Ia Id A d rho
global CLo CLa CDo CDa CMo CMa CRr
global CL_data CD_data CM_data CRr_rad CRr_AdvR CRr_data
global CMq CRp CNr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CoefUsed=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define non-aerodynamic parameters
m = 0.175; % Kg
g = 9.7935; % m/s^2
%A = 0.057; % m^2
A = .038; % disc golf (9" diameter)
d = 2*sqrt(A/pi); % diameter
rho = 1.23; % Kg/m^309o
Ia = 0.002352; % moment of inertia about the spinning axis
Id = 0.001219; % moment of inertia about the planar axis'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE THREE ESTIMATED COEFFICIENTS
%CMq= -0.005, CRp =-0.0055, CNr = 0.0000071 % short (three) flights
CMq= -1.44E-02; CRp =-1.25E-02; CNr = -3.41E-05; % long flight f2302

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

CRr_data = [-0.0172 -0.0192 -0.018 -0.0192 -0.018 -0.0172 -0.0172 -0.0168 -0.0188 -0.0164 -0.0136 -0.01 -0.0104 -0.0108 -0.0084 -0.008 -0.008 -0.006 -0.0048 -0.0064 -0.008 -0.003 ...
    -0.0112 -0.0132 -0.012 -0.0132 -0.012 -0.0112 -0.0112 -0.0108 -0.0128 -0.0104 -0.0096 -0.0068 -0.0072 -0.0076 -0.0052 -0.0048 -0.0048 -0.0028 -0.0032 -0.0048 -0.0064 -0.003 ...
    -0.0056 -0.0064 -0.0064 -0.0068 -0.0064 -0.0064 -0.0052 -0.0064 -0.0028 -0.0028 -0.004 -0.002 -0.004 -0.002 -0.0016 0 0 0 0 -0.002 -0.0048 -0.003 ...
    -0.0012 -0.0016 -0.0004 -0.0028 -0.0016 -0.0016 -0.0004 0.0004 0.0004 0.0008 0.0004 0.0008 0.0012 0.0008 0.002 0.0028 0.0032 0.0024 0.0028 0.0004 -0.0012 -0.003 ...
    -0.0012 -0.0012 -0.0016 -0.0016 -0.0012 -0.0004 0.0004 0.0008 0.0008 0.0016 0.0004 0.002 0.0004 0.0016 0.002 0.002 0.002 0.0012 0.0012 0 -0.0012 -0.003 ...
    -0.0012 -0.0012 -0.0004 -0.0008 -0.0008 -0.0008 0.0004 0.0004 0.0004 0.0008 0.0004 0.0008 -0.0004 0 0 0.0004 0 0 0.0004 -0.002 -0.0012 -0.003];
CRr_data = reshape(CRr_data,[22 6]);   

% CRr_deg=[-5 -4 -3 -2 -1 0 1 2 3 4
% 5 6 7 8 9 10 11 12 13 14 15 30 ]
CRr_rad = [-0.0873 -0.0698 -0.0524 -0.0349 -0.0175 0.0000 0.0175 0.0349 0.0524 0.0698 0.0873 0.1047 0.1222 0.1396 0.1571 0.1745 0.1920 0.2094 0.2269 0.2443 0.2618 0.5236];
CRr_AdvR= [2 1.04 0.69 0.35 0.17 0];

% % Define Simulation set initial conditions, enter your choosen values here:
% thetao = .192; % initial pitch angle
% speedo = 13.7; % magnitude, m/sec
% betao = .105; % flight path angle in radians between velocity vector and horizontal
% gd=50; % initial spin
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alphao = thetao - betao; % inital alpha
% vxo = speedo * cos(betao); % velocity component in the nx direction
% vzo = -(speedo * sin(betao)); % velocity component in the nz direction
% (note: nz is positive down)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0= vector of initial conditions
% x0= [ x y z vx vy vz phi theta phidot thetadot gd gamma]
x0=[-9.03E-01 -6.33E-01  -9.13E-01 1.34E+01 -4.11E-01 1.12E-03 -7.11E-02 2.11E-01 -1.49E+01 -1.48E+00 5.43E+01 5.03E+00];
x0(1:2)=0; %set origin to zero
% %x0(3) = -.9; % set height (negative z is positive... assholes)
% %x0(5)=0; % set y velocity to zero
% %x0(6) = x0(4).*-.1; % set z velocity
% %x0(4) = x0(4)*9; %scale throw speed
% %noseangle = 1; %degrees
% %x0(8) = noseangle.*pi/180; %pitch angle (nose up ness) in radians
 fps = 55; %throw speed fps
 x0(4) = fps/3.28; %intermediate throw speed
 x0(11) = x0(11)*1; % scale spin speed
% x0(7:10)=0; %throw flat with no off axis torque

CD_data = CD_data.*1;           %scale drag data to match a disc golf disc (lower than lids)

%offset = linspace(-.2,.5,6); %degrees
offset = [.3];
nR = length(offset);
nsteps = 300;
Paths = zeros(nR,3,nsteps+1); %[releaseangle x y]
ts = zeros([nR 1]);

for i=1:nR
i
    %change CRr
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CL_data=[ -0.1745 -0.2250 -10;
                -0.05236 0 -3;
                0 0.150 0;
                0.08727 0.4500 5;
                0.17453 0.7250 10;
                0.26180 0.9750 15;
                0.34907 1.2000 20;
                0.43633 1.4500 25;
                0.52360 1.6750 30];
    CL_data(:,2) = CL_data(:,2) + offset(i);
    
    % Enter values for tfinal and nsteps:
    tfinal = 6; % length of flight
    tspan=[0:tfinal/nsteps:tfinal];
    %options=[]
    options = odeset('AbsTol', 0.000001,'RelTol', 0.00000001,'Events',@hitZero);
    %'OutputFcn',@odeplot
    %options = odeset('AbsTol', 0.000001,'RelTol', 0.00000001);
    % Calls the ODE and integrate the frisbee equations of motions in the
    % subroutine, discfltEOM.m
    [t,x]=ode45(@discfltEOM,tspan,x0,options,CoefUsed);
    xx = x(:,1);
    yy = x(:,2);
    zz = x(:,3);
    Paths(i,1,1:length(xx)) = 3.28.*xx;
    Paths(i,2,1:length(yy)) = 3.28.*yy;
    Paths(i,3,1:length(zz)) = 3.28.*zz;
    ts(i)=length(t);
    
end

% Visualize Paths
%
figure
subplot(2,1,1)

for i=1:nR
    I = ts(i);
    hold on
    plot(squeeze(Paths(i,1,1:I)),squeeze(Paths(i,2,1:I)))
    axis equal
end
title('XY plane')

subplot(2,1,2)
for i=1:nR
    I = ts(i);
    hold on
    plot(squeeze(Paths(i,1,1:I)),-1.*squeeze(Paths(i,3,1:I)))
    axis equal
end
title('XZ plane')


% Event function to stop integration at z=0;

function [value, isterminal, direction] = hitZero(t,x,opArgs)
    value = x(3);
    isterminal = 1;
    direction = 0;
end

