% two mass connected by a spring and damper
clear;
% system definition
M=1; K=1; C=0.1;
A = [0 1 0 0; -K/M -C/M K/M C/M; 0 0 0 1; K/M C/M -K/M -C/M];
B = [0;1/M;0;0];
C = [0 0 1 0];
D =[0];
Sys =ss(A,B,C,D);
% Inverse System
By = (1/(C*A*A*B));
Ky = C*A*A*A/(C*A*A*B);
Ainv = A-(B*Ky); Binv = B*By; Cinv = -Ky; Dinv =By;
Sys_inv = ss(Ainv,Binv,Cinv,Dinv);
% defining the desired output and filtering it three times
tin = 0.5; tup = 1; tf = 5; delt = 0.001; ymax = 10;
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
Wf = 10 % choose a filter break frequency in Hz
num = [Wf*2*pi]; den = [1 (Wf*2*pi)]; % first order filter
[Af,Bf,Cf,Df] = tf2ss(num,den);
Sys_f = ss(Af,Bf,Cf,Df); % first order filter system
Sysf = Sys_f*Sys_f*Sys_f; % third order filter
[yd,xtemp]= lsim(Sysf,y,t); % filter the desired signal
figure(1); subplot(212), plot(t,yd);
xlabel('time'); ylabel('yd')
% numerical differentiation
y3d = diff(yd,3)/(delt^3); % the third derivative of yd
[mm,nn] = size(y3d);
time = 0:delt:(mm-1)*delt;
figure(2); clf; subplot(211), plot(time,y3d)
xlabel('time'); ylabel('y^{(3)}')
%The inverse input
[uff,xref] = lsim(Sys_inv,y3d,time);
figure(2); subplot(212), plot(time,uff)
xlabel('time'); ylabel('u_{ff}')
% Verify that this is correct...
[y_actual,x_actual]=lsim(Sys,uff,time);
figure(3); clf; plot(t,yd,'.g',time,y_actual,'r')
xlabel('time'); ylabel('y_{actual} in red')
legend('yd','y(actual)')
zoom on 