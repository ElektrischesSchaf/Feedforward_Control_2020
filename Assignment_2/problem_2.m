clc; clear;
r=2;

%% defining the desired output (and filtering it twice)
tin = 1; tup =3;  tf = 8; delt = 0.001; ymax = 10;

ramp = ymax/(tup-tin);
t1 = 0:delt:tin; 
t2=max(t1)+delt:delt:tup;
t3=max(t2)+delt:delt:tf;
y0 = zeros(size(t1));
y1 = ramp*(t2-max(t1));
y2 = max(y1)*ones(size(t3));
t = 0:delt:tf; y = [y0 y1 y2];

%% filtering the desired signal 
Wf = 1; % choose a filter break frequency in Hz
num = [Wf*2*pi]; den = [1 (Wf*2*pi)]; % first order filter
[Af,Bf,Cf,Df] = tf2ss(num,den);
Sys_f = ss(Af,Bf,Cf,Df); % first order filter system
Sysf = Sys_f*Sys_f*Sys_f*Sys_f*Sys_f; % fifth order filter
[yd,xtemp]= lsim(Sysf,y,t); % filter the desired signal
% figure(1);
% subplot(211), plot(t,y); 
% axis([0 8 -1 11]);
% xlabel('time(s)'); ylabel('y (unfiltered)')
% subplot(212), plot(t,yd); 
% xlabel('time(s)'); ylabel('yd (filtered)')
% axis([0 8 -1 11]);

%% yd_double_dot
for i=1:length(y)-1
    y_velocity(i)=y(i+1)-y(i);
end

for i=1:length(y_velocity)-1
    y_acc(i)=y_velocity(i+1)-y_velocity(i);
end

% https://stackoverflow.com/questions/26916066/variable-input-changing-in-time-ode45-matlab
eta1=0; eta2=0;
IC=[eta1, eta2];
time_span=[0 tf];
[time, eta] = ode45( @(time, eta) func(time, eta, y_acc), time_span, IC);

function eta_dot=func(time, eta, y_acc)
    M=1; m=1; L=1; B=1; g=9.8;
    eta_dot=zeros(2,1);
    
    eta_dot(1)=eta(2);
    eta_dot(2)=( -1/(m*L*L) )*(m*g*sin(eta(1))+B*eta(2)) + ((-1/m*L*L))*(m*L*cos(eta(1)))*y_acc(time);
end

