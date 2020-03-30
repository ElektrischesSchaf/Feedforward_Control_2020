clc; 
clear;

L=1; R=1; K=1; m=1; g=9.8; J=1; B=1;
%% origianl system state space
A=[
    0, 1, 0;
    0, -B/J, K/J;
    0, -K/L, -R/L;
    ];
B=[ 0 ; 0 ; 1/L];
C=[1, 0, 0];
D=0;
%% origianl system transfer function
[numerator denominator]=ss2tf(A,B,C,D);
r=3;
origianl_system_transfer_function=tf(numerator, denominator); % got relative degree=3
original_system_state_space=ss(A, B, C, D);

%% yd
% defining the desired output (and filtering it twice)
tin = 1; tup =3;  tf = 8; delt = 0.001; ymax = 10;

ramp = ymax/(tup-tin);
t1 = 0:delt:tin; 
t2=max(t1)+delt:delt:tup;
t3=max(t2)+delt:delt:tf;
y0 = zeros(size(t1));
y1 = ramp*(t2-max(t1));
y2 = max(y1)*ones(size(t3));
t = 0:delt:tf; y = [y0 y1 y2];

% filtering the desired signal 
Wf = 1; % choose a filter break frequency in Hz
num = [Wf*2*pi]; den = [1 (Wf*2*pi)]; % first order filter
[Af,Bf,Cf,Df] = tf2ss(num,den);
Sys_f = ss(Af,Bf,Cf,Df); % first order filter system
Sysf = Sys_f*Sys_f*Sys_f*Sys_f*Sys_f; % fifth order filter
[yd,xtemp]= lsim(Sysf,y,t); % filter the desired signal
figure(1);
subplot(211), plot(t,y); 
title('Desired Input y');
axis([0 8 -1 11]);
xlabel('time(s)'); ylabel('y (unfiltered)')

subplot(212), plot(t,yd); 
title('Desired Input yd');
xlabel('time(s)'); ylabel('yd (filtered)')
axis([0 8 -1 11]);

%% y_desired_double_dot
cc=[0]
for i=1:length(y)-1
    y_velocity(i)=yd(i+1)-yd(i);
end

y_velocity=[y_velocity cc];
for i=1:length(y_velocity)-1
    y_acc(i)=y_velocity(i+1)-y_velocity(i);
end
y_acc=[y_acc cc];

for i=1:length(y_acc)-1
    y_r(i)=y_acc(i+1)-y_acc(i);
end
y_r=[y_acc cc];

T=[ C ; C*A; C*A^(r-1) ];

Z_ref=[y ; y_velocity ; y_acc];

X_ref=inv(T)*Z_ref;

U_inv=( 1/ (C*A^(r-1)*B) )* ( y_r' -C*A^(r)* X_ref);

figure(2)
U_ff=U_inv(1,:); % TODO
subplot(211); plot(U_ff);
title('U_f_f');
y_angle=lsim(original_system_state_space, U_ff, t);
subplot(212); plot(y_angle);