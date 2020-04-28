clc;
clear;
LineWidth=3;
%% The relative degree of the system
% relative degree =  6-4 =2
r=2;
%% initialize the system
numerator=[11.88 4.977 539.6 129.9 5625]
denominator=[1 1.169 50.3 45.94 685.3 391.7 1952]
Origianl_System=tf(numerator, denominator);

%% Original System State Space 
% https://en.wikibooks.org/wiki/Control_Systems/Standard_Forms
Original_System_State_Space= canon(Origianl_System,'companion')
A=Original_System_State_Space.A';
B=Original_System_State_Space.C';
C=Original_System_State_Space.B';
D=Original_System_State_Space.D';

%% U_ff signal generation
%% 1-(e)
t=0:0.1:10;
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]); 
[tt,yy]=ode45(@findy,[0 10],[0 0],options);
yd=interp1(tt,yy(:,1),t);
yv=interp1(tt,yy(:,2),t);
figure(1)
hold on
plot(t, yd,'-')
plot(t, yv, '*')
title('Desired displacement and velocity profile');

%% Try reduce order method
T_t=[C; C*A];
T_b=[
   0 0 1 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   ];
T=[T_t; T_b];
Z_ref=[yd; yv];
T_inv=inv(T);
T_inv_L=T_inv(:, 1:2);
T_inv_R=T_inv(:, 3:6);
A_inv_red=T_b*A*T_inv_R-(T_b*B*C*(A^r)*T_inv_R)/(C*A^(r-1)*B);
B_inv_red=[ T_b*A*T_inv_L-(T_b*B*C*(A^r)*T_inv_L)/(C*A^(r-1)*B) ,  (T_b*B)/(C*A^(r-1)*B) ];
cc=[0]
for i=1:length(yv)-1
    ya(i)=yv(i+1)-yv(i);
end
ya=[ya cc];
Y_d=[yd; yv; ya]; % yd, yd_dot, yd_dot_dot
time_span=t;
IC=[0, 0, 0, 0];
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6]); % don't forget this !!
[t_out, eta_result]=ode45( @(t_in, eta) solve_eta(t_in, eta, A_inv_red, B_inv_red, Y_d), time_span, IC, options);
Z_ref=[ yd; yv; transpose(eta_result)];
X_ref=inv(T)*Z_ref;
U_inv_reduce_order=(1/(C*A^(r-1)*B))*(ya-C*A^(r)*X_ref);

figure(2);
plot(U_inv_reduce_order, 'LineWidth', LineWidth);
title('U_i_n_v by reduce order');

y_result=lsim(Original_System_State_Space, U_inv_reduce_order, t);
figure(3);
plot(y_result, 'LineWidth', LineWidth);
title('y result by reduce order');

for i=1:length(y_result)-1
    y_result_dot(i)=y_result(i+1)-y_result(i);
end
y_result_dot=[y_result_dot cc];

for i=1:length(y_result_dot)-1
    y_result_dot_dot(i)=y_result_dot(i+1)-y_result_dot(i);
end
y_result_dot_dot=[y_result_dot_dot cc];

for i=1:length(y_result)
      error_1(i)=y_result(i)-yd(i);
      error_2(i)=y_result_dot(i)-yv(i);
end

error=[error_1; error_2];

v=-[2 3];
for i = 1:length(U_inv_reduce_order)
    u_online(i)=U_inv_reduce_order(i) + inv(C*A^(r-1)*B)*v*error(:,i);
end

figure(4);
plot(u_online, 'LineWidth', LineWidth);
title('U_o_n_l_i_n_e');

%%  findy
function dy=findy(time,y)
    dy=zeros(2,1);
    dy(1)=y(2);
    if time<1
       dy(2)=1;
    elseif (time>=1 && time<3)
       dy(2)= -1;
    elseif (time>=3 && time<4)
        dy(2)=1;
    else
        dy(2)=0;
    end
end

%% Solving Eta (reduce order)
function eta_dot= solve_eta(t_in, eta, A_inv_red, B_inv_red, Y_d)
    time=[ 0:0.1:10];% must be consistant to real time span, not choose arbitrary value such as 1:1:101
    y_d_1=interp1(time,  Y_d(1,:), t_in);
    y_d_2=interp1(time,  Y_d(2,:), t_in);
    y_d_3=interp1(time,  Y_d(3,:), t_in);
    new_Y_d=[y_d_1; y_d_2; y_d_3];
    eta_dot=A_inv_red*eta+B_inv_red*new_Y_d;
end

