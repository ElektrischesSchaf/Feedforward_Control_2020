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
Original_System_State_Space= canon(Origianl_System,'companion')
A=Original_System_State_Space.A;
B=Original_System_State_Space.B;
C=Original_System_State_Space.C;
D=Original_System_State_Space.D;

%% inverse system State Space
A_inv=A-((B*C*A*A)/(C*A*B));
B_inv=(B)/(C*A*B);
C_inv=(-C*A*A)/(C*A*B);
D_inv=1/(C*A*B);

%% U_ff signal generation
% https://electrosome.com/signal-generation-in-matlab/
t=0:0.01:10;
x1=1*[t>=0];
x2=-2*[t>=1];
x3=2*[t>=3];
x4=-1*[t>=4];
yd=x1+x2+x3+x4;
% plot(t,input); % stem(n,x);
% title('Unit Step Signal - Cont');

%% Inverse System
Inverse_System_State_Space=ss(A_inv, B_inv, C_inv, D_inv);
%% 
U_ff=lsim(Inverse_System_State_Space, yd, t);
plot(U_ff, 'LineWidth', LineWidth);
hold on
pl=plot(yd, 'LineWidth', LineWidth*5);
pl.Color(4)=0.25;
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
% clf;
plot(y_acceleration,'r-' ,'LineWidth', LineWidth);


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
plot(y_acceleration_plus_5_percent, 'LineWidth', LineWidth);

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

% lgd=legend('voltage input (u_f_f)', 'y_d', 'y acceleration', 'y_+_5_%', 'y_-_5_%');
% lgd.FontSize=20;
% hold off

%% PID Feedback controller
Kp=0.75;
Ki=0.5;
Kd=0.5;
H=[1];
Controller=pid(Kp,Ki,Kd);
Origianl_System_minus_five_percent_PID=feedback(Controller*Origianl_System_minus_five_percent,H);
y_position_minus_5_percent_PID=lsim(Origianl_System_minus_five_percent_PID, U_ff, t);
y_position_minus_5_percent_PID=1e2.*y_position_minus_5_percent_PID;
    %% Differentiate y_position twice into y_acceleration
for i=1:length(t)-1
    y_velocity_minus_5_percent_PID(i)=y_position_minus_5_percent_PID(i+1)-y_position_minus_5_percent_PID(i);
end

for i=1:length(t)-2
    y_acceleration_minus_5_percent_PID(i)=y_velocity_minus_5_percent_PID(i+1)-y_velocity_minus_5_percent_PID(i);
end
plot(y_acceleration_minus_5_percent_PID, 'LineWidth', LineWidth);
set(gcf,'Position',[50 50 750 750])
lgd=legend('voltage input (u_f_f)', 'y_d', 'y acceleration', 'y_+_5_%', 'y_-_5_%', '-5% with PID');
lgd.FontSize=20;
