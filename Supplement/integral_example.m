% This example shows how to find the integral with command ode45
% Szu-Chi Tien 2008
clear all;
close all;


%set second derivative of desired trjectory (i.e. y-double-dot)
delt=0.001; 
t1=0:delt:1;                  % set four sectors
t2=1+delt:delt:3;
t3=3+delt:delt:4;
t4=4+delt:delt:10;
y2d_1=ones(size(t1));
y2d_2=-1*ones(size(t2));
y2d_3=ones(size(t3));
y2d_4=zeros(size(t4));
t=0:delt:10; 
y2d=[y2d_1 y2d_2 y2d_3 y2d_4];  % this is the second derivative of desired trajectory

% get the desired trjectory
yd0=0; % initial value of yd
ydot0=0; % initial value of y-dot
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
[time,yd_temp]=ode45('findy',[0 10],[yd0 ydot0],options);
yd=zeros(size(t));
for i=1:length(t),
    yd(i)=interp1(time,yd_temp(:,1),delt*(i-1));
end

figure(1)
plot(t, y2d);
xlabel('time(ms)');
ylabel('y^{(2)}');
axis([0 10 -1.5 1.5]);


figure(2)
plot(t,yd,'b-');
xlabel('time(ms)');
ylabel('y');
