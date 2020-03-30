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

origianl_system_transfer_function=tf(numerator, denominator); % got relative degree=3
Original_System_State_Space=ss(A, B, C, D);
% Original_System_State_Space = canon(origianl_system_transfer_function,'companion')

A_o=Original_System_State_Space.A';
B_o=Original_System_State_Space.C';
C_o=Original_System_State_Space.B';
D_o=Original_System_State_Space.D';
origianl_system_transfer_function_2=ss2tf(A_o, B_o, C_o, D_o);

%% inverser system state space
r=3;
A_inv=A_o-((B_o*C_o*A_o^(r))/(C_o*A_o^(r-1)*B_o));
B_inv=(B_o)/(C_o*A_o^(r-1)*B_o);
C_inv=-(C_o*A_o*A_o)/(C_o*A_o^(r-1)*B_o);
D_inv=1/(C_o*A_o^(r-1)*B_o);


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


%% calculate inveres response
inverse_system_state_space=ss(A_inv, B_inv, C_inv, D_inv);
inverse_system_transfer_function=ss2tf(A_inv, B_inv, C_inv, D_inv);

U_ff=lsim(inverse_system_state_space, y, t);

my_y=lsim(Original_System_State_Space, U_ff, t);
figure(2);

subplot(211); plot(my_y);
title('My y');

subplot(212); plot(U_ff);
title('U_f_f');