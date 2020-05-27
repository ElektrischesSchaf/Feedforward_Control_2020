clc; clear;
r=2;
M=1; m=1; L=1; B=1; g=1;
%% defining the desired output (and filtering it twice)
tin = 7 ; tup = 12 ; tf = 22 ; delt = 0.001 ; ymax = 10;
ramp = ymax/(tup-tin);
t1 = 0:delt:tin; 
t2=max(t1)+delt:delt:tup;
t3=max(t2)+delt:delt:tf;
y0 = zeros(size(t1));
y1 = ramp*(t2-max(t1));
y2 = max(y1)*ones(size(t3));
t = 0:delt:tf; y = [y0 y1 y2];
figure(1); clf; subplot(211), plot(t,y); 
axis([0,tf,-0.1,1.2*ymax])
xlabel('time'); ylabel('y unfiltered')
% filtering the desired signal 
Wf = 1; % choose a filter break frequency in Hz
num = [Wf*2*pi]; den = [1 (Wf*2*pi)]; % first order filter
[Af,Bf,Cf,Df] = tf2ss(num,den);
Sys_f = ss(Af,Bf,Cf,Df); % first order filter system
Sysf = Sys_f*Sys_f*Sys_f*Sys_f*Sys_f; % fifth order filter
[yd,xtemp]= lsim(Sysf,y,t); % filter the desired signal
figure(1);  subplot(212), plot(t,yd); 
xlabel('time'); ylabel('yd')
% numerical differentiation
y1d = diff(yd,1)/(delt^1); % the first derivative of yd
y1d = [y1d;0];
y2d = diff(yd,2)/(delt^2); % the second derivative of yd
y2d = [y2d;0;0];
[mm,nn] = size(y2d);
time = 0:delt:(mm-1)*delt;
figure(1); clf; subplot(311), plot(time,yd)
xlabel('time'); ylabel('yd')
subplot(312), plot(time,y1d)
xlabel('time'); ylabel('yd^1')
subplot(313), plot(time,y2d)
xlabel('time'); ylabel('yd^2')


%% y_desired_double_dot
for i=1:length(y)-1
    y_velocity(i)=yd(i+1)-yd(i);
end
cc=[0]
y_velocity=[y_velocity cc];
for i=1:length(y_velocity)-1
    y_acc(i)=y_velocity(i+1)-y_velocity(i);
end
y_acc=[y_acc cc];

%% Initiate the iterative procedure 
% https://stackoverflow.com/questions/26916066/variable-input-changing-in-time-ode45-matlab
eta = zeros(2,length(time))';
theta = eta(:,1);

% time_span=[0 tf]; % If with this, the response will be downsampled.
time_span=t;
% options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]); 
% options2 = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6]); 
% [time, eta] = ode45( @(time, eta) func(time, eta, y_acc, M, m, L, B, g), time_span, IC, options);
Aold=[0 1;
   g/L -B/(m*L^2)];
[V,D]=eig(Aold);
TransM=inv(V);
check_1=TransM*Aold*inv(TransM);  % check internal dynamics decoupling


%% iterative procedure 
figure(3); clf
theta_store = [];
theta1d_store = [];
nsteps = 30;
for jj = 1:1:nsteps
    p=(g/L).*(sin(theta)-theta)+(1/L).*cos(theta).*y2d;
    [zs,zstemp] = lsim(-1.618, 0.8507, 1, 0, p, t);
    [zu,zutemp] = lsim(-0.618, -0.5257, 1, 0, flipud(p), t);

    %%%% compute the new internal dynamics
    eta = (inv(TransM)*[zs zu]')';
    error = norm(theta -eta(:,1))
    theta = eta(:,1);
    theta_1d = eta(:,2);
    theta_store = [theta_store theta];
    theta1d_store = [theta1d_store theta_1d];
    %%%% plot theta and internal states every iteration 
    subplot(311), plot(time,zs); hold on
    xlabel('time'); ylabel('z_s')
    subplot(312), plot(time,zu)
    xlabel('time'); ylabel('z_u'); hold on
    subplot(313), plot(time,theta)
    xlabel('time'); ylabel('\theta')
    drawnow
    zoom on; hold on
end
figure(3); hold off

figure(4); clf
subplot(211);
plot(time,theta_store,time,theta,'r')
xlabel('time'); ylabel('iterations of \theta')
subplot(212);
plot(time,theta1d_store,'b')
xlabel('time'); ylabel('iterations of \theta1d')

%{
M=1; m=1; L=1; B=1; g=9.8;
ttt= [0: delt : tf];
y_desired_dot_dot=interp1( ttt, y_acc, time );

alpha_2_dot=(-1/m*L^2)*( -m*g*L.*sin(eta(:,1))+B.*eta(:,2) - m*L.*cos(eta(:,1)).*y_desired_dot_dot );
U_inv =  (M+m).*y_desired_dot_dot  -  m*L.*cos(eta(:,1)).*(  alpha_2_dot )  +  m*L*sin(eta(:,1)).*eta(:,2).^2  ;


% yee=0:0.2:8;
% U_inv=interp1( 1:41, U_inv, 0:delt:8);
figure(2);
plot(U_inv);
title('U_f_f');
U=[U_inv ,  t'];

figure(3);
[t, x_n]=ode45(@(t, x_n) cart(t, x_n, U, M, m, L, M, g), time_span, [0 0 0 0], options2);
plot(x_n(:,1));
title('Origianl system output after applied U_f_f');



function eta_dot=func(time, eta, y_acc, M, m, L, B, g)
    
    ttt= [ 0: 0.001 : 22];
    y_desired_dot_dot=interp1( ttt, y_acc, time );
    
    eta_dot=zeros(2,1);
    
    eta_dot(1)=eta(2);
    eta_dot(2)=( -1/(m*L*L) )*(m*g*sin(eta(1))+B*eta(2)) + ((-1/m*L*L))*(m*L*cos(eta(1)))*y_desired_dot_dot;
end


% x_n(1)= Xc, x_n(2)= Xc_dot, x_n(3)= theta, x_n(4)= theta_dot
function x_diff = cart(t, x_n , U, M, m, L, B, g)
    
    F=U(:,1);
    time=U(:,2);    
    u = interp1(time, F, t);
    
    x_diff=zeros(2,1);
    x_diff(1)=x_n(2);
    x_diff(3)=x_n(4);
    A_s=[m*L*L  , -m*L.*cos(x_n(3));
        -m*L.*cos(x_n(3))   ,  M+m]
    A_s=A_s.*(1/(m*M*L*L+m*m*L*L.*sin(x_n(3)).*sin(x_n(3)) ) )
    B_s=[m*L.*sin(x_n(4)).*x_n(4)^2+u;
        -m*g*L.*sin(x_n(3))-B.*x_n(4)];
    C_s=A_s*B_s;
    x_diff(2)=C_s(1);
    x_diff(4)=C_s(2);

end

%}
