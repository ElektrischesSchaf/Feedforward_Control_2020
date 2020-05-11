clc;
clear;

A=[
    0 1 0;
    -1 -1 0;
    0 0 -1;
    ];
B=[
    1 0;
    0 1;
    1 0
    ];
C=[
    1 0 0;
    0 0 1;
    ];
D=0;
new_C1= [0 1 0;[0 -1 -1]*A]
new_D1=[1 0;[0 -1 -1]*B]

A_inv=A-B*inv(new_D1)*new_C1;
B_inv=B*inv(new_D1);
C_inv=-inv(new_D1)*new_C1;
D_inv=new_D1;

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

Y_d=[y_velocity; y_velocity.*0 ];
Silver_Inverse=ss(A_inv,B_inv,C_inv,D_inv);
U_ff=lsim(Silver_Inverse, Y_d, t);
figure(2);
plot(U_ff);
title('U_f_f by Silver method');

Original_System=ss(A,B,C,D)
y=lsim(Original_System, U_ff, t);
figure(3);
plot(y);
title('y by Silver method');