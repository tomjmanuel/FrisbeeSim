% Crowther, W. J., and J. R. Potts. "Simulation of a 
% spinstabilised sports disc. Sports Engineering 2007
% 10/5/19

% launch conditions (analysis of typical tragectory):

% displacement in earth axis: [m]
Xl=  -.9; Yl = -.63; Zl = -.91;

% body axis velocity components [m/s]
ul = 15; vl = 0; wl = 0;

% roll pitch and yawl euler angles at [rad]
phiL = 0; thetaL = 15; psiL = 0;

% body axis roll, pitch and yaw rates [rad/s]
%rl = 52.85; pl = -26.25; ql = -5.19;
rl = 0; pl = 0; ql = 0;

% spin rate [rev/s]
spinrate = 5;
omega = 2*pi*spinrate;

% mass
m = .175; %[kg]

% diameter
c = .275; %[m]
S = pi*(c/2)^2; %Area [m^2]

% get Force coefficients based on angle of attack
% taken from wind tunnel data
a = thetaL; %angle of attack
Clift = .13+3.09*a;
Cdrag = .085+3.3*(a+.052)^2;
Cm = -.01+.057*a; %pitching moment coeff

% First, solve F4
% F4 = qinf * S * Cf
% where Cf = [-Cdrag Cside -Clift]
% and qinf  1/2 (rho) * mag (X1dot)
















