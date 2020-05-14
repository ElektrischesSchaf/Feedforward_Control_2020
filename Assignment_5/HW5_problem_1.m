clc; clear;
[g_1_num, g_1_de] = zp2tf([-10 -20 -2]', [-3 -4 -5 -6 -7], 1);
[g_2_num, g_2_de] = zp2tf([-2]', [-3 -4 -5 -6 -7], 200);
[A_1, B_1, C_1, D_1] = tf2ss(g_1_num, g_1_de);
sys_1=ss(A_1,B_1,C_1,D_1);
[A_2, B_2, C_2, D_2] = tf2ss(g_2_num, g_2_de);
sys_2=ss(A_2,B_2,C_2,D_2);
r_1= 2; r_2=4;

%% defining the desired output and filtering it two times
tin = 1; tup = 3;  tdown=5; tback=7; tf = 10; delt = 0.001; ymax = 1;
ramp_1 = ymax/(tup-tin);
ramp_2 = -ymax/(tback-tdown);
t1 = 0:delt:tin; 
t2=max(t1)+delt:delt:tup;
t3=max(t2)+delt:delt:tdown;
t4=max(t3)+delt:delt:tback;
t5=max(t4)+delt:delt:tf;

y0 = zeros(size(t1));
y1 = ramp_1*(t2-max(t1));
y2 = max(y1)*ones(size(t3-max(t2)));
y3 = max(y2)+ ramp_2*(t4-max(t3));
y4 = zeros(size(t5-max(t4)));

t = 0:delt:tf;
y = [y0 y1 y2 y3 y4];
figure(1); clf; subplot(211), plot(t,y); 
axis([0, tf, -0.1, 1.2*ymax])
xlabel('time'); ylabel('y unfiltered')
% filtering the desired signal 
Wf = 1 % choose a filter break frequency in Hz
num = [Wf*2*pi]; den = [1 (Wf*2*pi)]; % first order filter
[Af,Bf,Cf,Df] = tf2ss(num,den);
Sys_f = ss(Af,Bf,Cf,Df); % first order filter system
Sysf = Sys_f*Sys_f*Sys_f*Sys_f; % second order filter
[yd,xtemp]= lsim(Sysf,y,t); % filter the desired signal
figure(1);  subplot(212), plot(t,yd); 
xlabel('time'); ylabel('yd')
%% numerical differentiation
y1d = diff(yd,1)/(delt^1); % the first derivative of yd
y1d = [y1d;0];
y2d = diff(yd,2)/(delt^2); % the second derivative of yd
y2d = [y2d;0;0];
y3d = diff(yd,3)/(delt^3); % the third derivative of yd
y3d = [y3d;0;0;0];
y4d = diff(yd,4)/(delt^4); % the fourth derivative of yd
y4d = [y4d;0;0;0;0];
[mm,nn] = size(y2d);
time = 0:delt:(mm-1)*delt;
figure(1); clf; 
subplot(511), plot(time,yd);xlabel('time'); ylabel('yd')
subplot(512), plot(time,y1d);xlabel('time'); ylabel('yd^1')
subplot(513), plot(time,y2d);xlabel('time'); ylabel('yd^2')
subplot(514), plot(time,y3d);xlabel('time'); ylabel('yd^3')
subplot(515), plot(time,y4d);xlabel('time'); ylabel('yd^4')

%% G1
T_t=[C_1; C_1*A_1];
T_b=[
   0 0 1 0 0 ;
   0 0 0 1 0 ;
   0 0 0 0 1 ;
   ];
T=[T_t; T_b];
T_inv=inv(T);
T_inv_L=T_inv(:, 1:2);
T_inv_R=T_inv(:, 3:5);
A_inv_red = T_b*A_1*T_inv_R-(T_b*B_1*C_1*(A_1^r_1)*T_inv_R)/(C_1*A_1^(r_1-1)*B_1);
B_inv_red = [ T_b*A_1*T_inv_L-(T_b*B_1*C_1*(A_1^r_1)*T_inv_L)/(C_1*A_1^(r_1-1)*B_1) ,  (T_b*B_1)/(C_1*A_1^(r_1-1)*B_1) ];
Y_d=[yd y1d y2d];
time_span=t;
IC=[0, 0, 0];
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 ]); % don't forget this !!
[t_out, eta_result]=ode45( @(t_in, eta) solve_eta_1(t_in, eta, A_inv_red, B_inv_red, Y_d), time_span, IC, options);
Z_ref=[ y1d y2d eta_result];
X_ref=inv(T)*Z_ref';
U_inv_reduce_order_1=(1/(C_1*A_1^(r_1-1)*B_1))*(y2d-C_1*A_1^(r_1)*X_ref);
y_result_1=lsim(sys_1, U_inv_reduce_order_1(:,1), t);

%% G2 relative degree=4
T_t=[C_1; C_1*A_1; C_1*A_1*A_1 ;C_1*A_1*A_1*A_1];
T_b=[
   0 0 1 0 0 
   ];
T=[T_t; T_b];
inv(T);
T_inv_L=T_inv(:, 1:4);
T_inv_R=T_inv(:, 5);
A_inv_red = T_b*A_1*T_inv_R-(T_b*B_1*C_1*(A_1^r_1)*T_inv_R)/(C_1*A_1^(r_1-1)*B_1);
B_inv_red = [ T_b*A_1*T_inv_L-(T_b*B_1*C_1*(A_1^r_1)*T_inv_L)/(C_1*A_1^(r_1-1)*B_1) ,  (T_b*B_1)/(C_1*A_1^(r_1-1)*B_1) ];
Y_d=[yd y1d y2d y3d y4d];
time_span=t;
IC=[0];
options = odeset('RelTol',1e-6,'AbsTol',[1e-6]); % don't forget this !!
[t_out, eta_result]=ode45( @(t_in, eta) solve_eta_2(t_in, eta, A_inv_red, B_inv_red, Y_d), time_span, IC, options);
Z_ref=[ yd y1d y2d y3d eta_result];
X_ref=inv(T)*Z_ref';
U_inv_reduce_order_2=(1/(C_2*A_2^(r_2-1)*B_2))*(y4d-C_2*A_2^(r_2)*X_ref);
y_result_2=lsim(sys_2, U_inv_reduce_order_2(:,1), t);

%% plot the result
figure(3); clf; 
subplot(311), plot(time,yd);xlabel('time'); ylabel('yd');
subplot(312), plot(time,y_result_1);xlabel('time'); ylabel('y with u_{ff}');
subplot(313), plot(time,y_result_2);xlabel('time'); ylabel('y with u_{ff}');

figure(4); clf;
plot(time,y_result_1);
xlabel('time'); ylabel('y with u_{ff} in sys 1');

figure(5); clf;
plot(time,y_result_2);
xlabel('time'); ylabel('y with u_{ff} in sys 2');
%% Solving Eta (reduce order)
function eta_dot= solve_eta_1(t_in, eta, A_inv_red, B_inv_red, Y_d)
    time=[0:0.001:10];% must be consistant to real time span, not choose arbitrary value such as 1:1:101
    y_d_1=interp1(time,  Y_d(:,1), t_in);
    y_d_2=interp1(time,  Y_d(:,2), t_in);
    y_d_3=interp1(time,  Y_d(:,3), t_in);
    new_Y_d=[y_d_1; y_d_2; y_d_3];
    eta_dot=A_inv_red*eta+B_inv_red*new_Y_d;
end


function eta_dot= solve_eta_2(t_in, eta, A_inv_red, B_inv_red, Y_d)
    time=[0:0.001:10];% must be consistant to real time span, not choose arbitrary value such as 1:1:101
    y_d_1=interp1(time,  Y_d(:,1), t_in);
    y_d_2=interp1(time,  Y_d(:,2), t_in);
    y_d_3=interp1(time,  Y_d(:,3), t_in);
    y_d_4=interp1(time,  Y_d(:,4), t_in);
    y_d_5=interp1(time,  Y_d(:,5), t_in);
    new_Y_d=[y_d_1; y_d_2; y_d_3; y_d_4; y_d_5];
    eta_dot=A_inv_red*eta+B_inv_red*new_Y_d;
end
