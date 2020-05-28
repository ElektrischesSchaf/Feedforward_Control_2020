clc; clear;
A=[
    0 1 0 0;
    1 0 1 0;
    0 0 0 1;
    0 0 0 0;
    ];
B=[1 0 ; 0 0 ; 0 1 ; 1 0];
C=[1 0 0 0;
    0 0 0 1];
C_1=C(1,:); C_2=C(2,:);
D=[0];
sys=ss(A,B,C,D);
% [num1 den1] = ss2tf(A,B,C,D,1);
% [num1 den1] = ss2tf(A,B,C,D,1);  % iu = 1
% [num2 den2] = ss2tf(A,B,C,D,2);  % iu = 2
tf(sys);
r_1=3; r_2=1;
C_2*A^1*B;
Ay=[ C_1*A^r_1 ; C_2*A^r_2 ];
By=[ C_1*A^(r_1-1)*B ; C_2*A^(r_2-1)*B];
inv(By)
T_t=[ C_1; C_1*A; C_2 ];
T_b=[0 0 1 0];
T=[T_t ; T_b];
T_inv=inv(T);
T_inv_L=T_inv(:, 1:3);
T_inv_R=T_inv(:, 4);
A_inv_red = T_b*A*T_inv_R - (T_b*B *inv(By)* Ay*T_inv_R);
B_inv_red = [ T_b*A*T_inv_L - (T_b*B* inv(By)* Ay *T_inv_L) ,  T_b*B*inv(By) ];
C_inv_red = -inv(By)*Ay*T_inv_R;
D_inv_red = [-inv(By)*Ay*T_inv_L , inv(By) ];

%% defining the desired output (and filtering it twice)
tin = 1 ; tup = 3 ; tf = 8 ; delt = 0.001 ; ymax = 10;
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
y3d = diff(yd,3)/(delt^3); % the second derivative of yd
y3d = [y3d;0;0;0];
[mm,nn] = size(y2d);
time = 0:delt:(mm-1)*delt;
%{
figure(1); clf; subplot(311), plot(time,yd)
xlabel('time'); ylabel('yd')
subplot(312), plot(time,y1d)
xlabel('time'); ylabel('yd^1')
subplot(313), plot(time,y2d)
xlabel('time'); ylabel('yd^2')
%}

Y_d=[ yd, y1d , yd, y3d ,y1d ]
u_inv= D_inv_red*Y_d';
input_1=u_inv(1,:);
input_2=u_inv(2,:);
figure(2); clf; 

subplot(211), plot(time, input_1);
title('U_{ff}');
subplot(212), plot(time, input_2);

[my_output, x_ref]= lsim( sys, u_inv, time);

figure(3); clf;
plot(time, my_output);
title('desired outpus');