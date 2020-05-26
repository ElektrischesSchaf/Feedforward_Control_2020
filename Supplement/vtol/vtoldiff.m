% this is the nonlinear model of the VTOL aircraft 

function xdiff = fun(t,x,flag,param,U)
%t  % print out simulation time to check progress of simulation

% unravel the parameteres passed into ode 
time = U(:,1); U1 = U(:,2); U2 = U(:,3); 
Thetad = U(:,4); Theta_1d = U(:,5);
lambda = param(1); epsilon = param(2); g = param(3); gamma = param(4);

u1 = interp1(time,U1,t); 
u2 = interp1(time,U2,t);
thetad = interp1(time,Thetad,t);
theta_1d = interp1(time,Theta_1d,t);

theta = x(5);

% The feedback controller PD 
Kp=0.01;
Kd=0.01;
u2_fb = -Kp*(x(5)-thetad) -Kd*(x(6)-theta_1d);
u2 = u2 + u2_fb;

xdiff = zeros(6,1);
xdiff(1) = x(2); 
xdiff(2) = -u1*sin(theta) +epsilon*u2*cos(theta); 
xdiff(3) = x(4); 
xdiff(4) = u1*cos(theta) +epsilon*u2*sin(theta) -g;
xdiff(5) = x(6);
xdiff(6) = lambda*u2;

return