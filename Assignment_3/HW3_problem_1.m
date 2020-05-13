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
% Z_ref=[yd; yv]; % wrong
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
Y_d=[yd; yv; ya];
time_span=t;
IC=[0, 0, 0, 0];
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6]); % don't forget this !!
[t_out, eta_result]=ode45( @(t_in, eta) solve_eta(t_in, eta, A_inv_red, B_inv_red, Y_d), time_span, IC, options);
Z_ref=[ yd; yv; transpose(eta_result)];
X_ref=inv(T)*Z_ref;
U_inv_reduce_order=(1/(C*A^(r-1)*B))*(ya-C*A^(r)*X_ref);



%% inverse system State Space
A_inv=A-((B*C*A*A)/(C*A*B));
B_inv=(B)/(C*A*B);
C_inv=(-C*A*A)/(C*A*B);
D_inv=1/(C*A*B);

%% Inverse System
Inverse_System_State_Space=ss(A_inv, B_inv, C_inv, D_inv);
Inverse_System=ss2tf(A_inv, B_inv, C_inv, D_inv);

%% 
t=0:0.01:10;
y_desired=[1*ones(1,100) -1*ones(1,200) 1*ones(1,100) 0*ones(1,601)];
[U_ff, t_inverse_system, X_ref]=lsim(Inverse_System_State_Space, y_desired, t);

figure(2);
plot(U_ff, 'LineWidth', LineWidth);
xlabel('time(ms)');
ylabel('Voltage');
title('The Inverse System Input');
figure(3)
y_result=lsim(Original_System_State_Space, U_ff, t);
plot(y_result, 'LineWidth', LineWidth);
title('The y response of the original system +5%');
xlabel('time(ms)');
ylabel('y (m)');

%% PID Feedback controller
Kp=3;
Ki=5;
Kd=0.5;
H=[1];
Controller=pid(Kp,Ki,Kd);
% Origianl_System_minus_five_percent_PID=U_ff+feedback(Controller*Origianl_System_minus_five_percent, H); % inverse and feedback input
PID=feedback(Controller*Origianl_System, H); % only feedback input
feedback_y=lsim(PID, U_ff, t);
figure(4);
plot(feedback_y, 'LineWidth', LineWidth);
title('Feedback control loop response');

%% Close loop
K=0.01;
A_cl=(A-B*K);
B_cl=B;
C_cl=C;
D_cl=D;

input_close_loop=U_ff+K*X_ref;

system_close_loop=ss( A_cl, B_cl, C_cl, D_cl );
y_close_loop_position=lsim(system_close_loop, input_close_loop(:,1), t);

figure(5);
plot(y_close_loop_position, 'LineWidth', LineWidth);
title('Feedforward + feedback control loop response');

figure(6);
plot(U_inv_reduce_order, 'LineWidth', LineWidth);
title('U_i_n_v by reduce order');

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

