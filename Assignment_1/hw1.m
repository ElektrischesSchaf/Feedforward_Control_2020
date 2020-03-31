clc;
clear;
LineWidth=3;
%% initialize the system
numerator=[11.88 4.977 539.6 129.9 5625]
denominator=[1 1.169 50.3 45.94 685.3 391.7 1952]
Origianl_System=tf(numerator, denominator);

%% Bode Plot of the system
% bode(sys); % change frequency unit to KHz in "View" -> "Property Editor"

%% Step response of the system
% plot(step(sys));

%% The relative degree of the system
% relative degree =  6-4 =2

%% Original System State Space 
% https://en.wikibooks.org/wiki/Control_Systems/Standard_Forms
Original_System_State_Space= canon(Origianl_System,'companion')
A=Original_System_State_Space.A';
B=Original_System_State_Space.C';
C=Original_System_State_Space.B';
D=Original_System_State_Space.D';

%% inverse system State Space
A_inv=A-((B*C*A*A)/(C*A*B));
B_inv=(B)/(C*A*B);
C_inv=(-C*A*A)/(C*A*B);
D_inv=1/(C*A*B);

%% U_ff signal generation
% https://electrosome.com/signal-generation-in-matlab/
start_time=0;
delta=0.01;
end_time=10
t=start_time:delta:end_time;
x1=1*[t>=0];
x2=-2*[t>=1];
x3=2*[t>=3];
x4=-1*[t>=4];
yd=x1+x2+x3+x4;
figure(1);
subplot(2,1,1);
plot(t,yd, 'LineWidth',LineWidth); % stem(n,x);
xlabel('time(ms)');
ylabel('V/ms^2');
title('The desired acceleration profile');


%% Inverse System
Inverse_System_State_Space=ss(A_inv, B_inv, C_inv, D_inv);
Inverse_System=ss2tf(A_inv, B_inv, C_inv, D_inv);
%% 
[U_ff, t_inverse_system, X_ref]=lsim(Inverse_System_State_Space, yd, t);
subplot(2,1,2)
plot(U_ff, 'LineWidth', LineWidth);
xlabel('time(ms)');
xlim([start_time end_time/delta]);
ylabel('Voltage');
title('The Inverse System Input');

% pl=plot(yd, 'LineWidth', LineWidth*5);
% pl.Color(4)=0.25;
%%
y_position=lsim(Original_System_State_Space, U_ff, t);
% plot(y_position);
y_position=1e4*y_position;

%% Differentiate y_position twice into y_acceleration
for i=1:length(t)-1
    y_velocity(i)=y_position(i+1)-y_position(i);
end

for i=1:length(t)-2
    y_acceleration(i)=y_velocity(i+1)-y_velocity(i);
end

figure(2);
plot(yd, 'LineWidth',LineWidth/2); % stem(n,x);
hold on;
plot(y_acceleration,'r-' ,'LineWidth', LineWidth/2);
xlabel('time(ms)');
xlim([start_time end_time/delta]);
ylabel('V/ms^2');
title('y_d and y_acceleration');
lgd=legend('y_d', 'y_acceleration');
lgd.FontSize=20;



%% Numerator  constant 5625 has changed by 5% 
numerator=[11.88 4.977 539.6 129.9 5625*1.05];
denominator=[1 1.169 50.3 45.94 685.3 391.7 1952];
Origianl_System_plus_five_percent=tf(numerator, denominator);
y_position_plus_5_percent=lsim(Origianl_System_plus_five_percent, U_ff, t);
y_position_plus_5_percent=1e4.*y_position_plus_5_percent;
    %% Differentiate y_position twice into y_acceleration
for i=1:length(t)-1
    y_velocity_plus_5_percent(i)=y_position_plus_5_percent(i+1)-y_position_plus_5_percent(i);
end

for i=1:length(t)-2
    y_acceleration_plus_5_percent(i)=y_velocity_plus_5_percent(i+1)-y_velocity_plus_5_percent(i);
end
figure(3);
plot(y_acceleration_plus_5_percent, 'LineWidth', LineWidth);
hold on;

%% Numerator  constant 5625 has changed by -5% 
numerator=[11.88 4.977 539.6 129.9 5625*0.95];
denominator=[1 1.169 50.3 45.94 685.3 391.7 1952];
Origianl_System_minus_five_percent=tf(numerator, denominator);
y_position_minus_5_percent=lsim(Origianl_System_minus_five_percent, U_ff, t);
y_position_minus_5_percent=1e4.*y_position_minus_5_percent;
    %% Differentiate y_position twice into y_acceleration
for i=1:length(t)-1
    y_velocity_minus_5_percent(i)=y_position_minus_5_percent(i+1)-y_position_minus_5_percent(i);
end

for i=1:length(t)-2
    y_acceleration_minus_5_percent(i)=y_velocity_minus_5_percent(i+1)-y_velocity_minus_5_percent(i);
end

plot(y_acceleration_minus_5_percent, 'LineWidth', LineWidth);
xlabel('time(ms)');
xlim([start_time end_time/delta]);
ylabel('V/ms^2');
title('+5% and -5% comparison');
lgd=legend('+5%', '-5%');
lgd.FontSize=20;

close all;

%% PID Feedback controller
Kp=30;
Ki=50;
Kd=0.5;
H=[1];
Controller=pid(Kp,Ki,Kd);
% Origianl_System_minus_five_percent_PID=U_ff+feedback(Controller*Origianl_System_minus_five_percent, H); % inverse and feedback input
Origianl_System_minus_five_percent_PID=feedback(Controller*Origianl_System_minus_five_percent, H); % only feedback input
y_position_minus_5_percent_PID=lsim(Origianl_System_minus_five_percent_PID, U_ff, t);

y_position_minus_5_percent_PID=1e2.*y_position_minus_5_percent_PID;
    %% Differentiate y_position twice into y_acceleration
for i=1:length(t)-1
    y_velocity_minus_5_percent_PID(i)=y_position_minus_5_percent_PID(i+1)-y_position_minus_5_percent_PID(i);
end

for i=1:length(t)-2
    y_acceleration_minus_5_percent_PID(i)=y_velocity_minus_5_percent_PID(i+1)-y_velocity_minus_5_percent_PID(i);
end


figure(4);
plot(yd, 'LineWidth',LineWidth/2); % stem(n,x);
hold on;
plot(y_acceleration_minus_5_percent_PID, 'LineWidth', LineWidth);
lgd=legend('y_d', 'PID feedback response');
lgd.FontSize=10;

%% Close loop
K=0.0;
A_cl=(A-B*K);
B_cl=B;
C_cl=C;
D_cl=D;

input_close_loop=U_ff+K*X_ref;

system_close_loop=ss( A_cl, B_cl, C_cl, D_cl );
y_close_loop_position=lsim(system_close_loop, input_close_loop(:,1), t);



for i=1:length(t)-1
    y_close_loop_velocity(i)=y_close_loop_position(i+1)-y_close_loop_position(i);
end

for i=1:length(t)-2
    y_close_loop_acceleration(i)=y_close_loop_velocity(i+1)-y_close_loop_velocity(i);
end

figure(5);
plot(y_close_loop_acceleration, 'LineWidth', LineWidth);
lgd=legend('ff and fd control response');
lgd.FontSize=10;