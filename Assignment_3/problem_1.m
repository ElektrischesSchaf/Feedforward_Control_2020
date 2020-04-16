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

%% inverse system State Space
A_inv=A-((B*C*A*A)/(C*A*B));
B_inv=(B)/(C*A*B);
C_inv=(-C*A*A)/(C*A*B);
D_inv=1/(C*A*B);

%% Inverse System
Inverse_System_State_Space=ss(A_inv, B_inv, C_inv, D_inv);
Inverse_System=ss2tf(A_inv, B_inv, C_inv, D_inv);

%% 
[U_ff, t_inverse_system, X_ref]=lsim(Inverse_System_State_Space, yd, t);

figure(2);
plot(U_ff, 'LineWidth', LineWidth);
xlabel('time(ms)');
ylabel('Voltage');
title('The Inverse System Input');
hold on;
y_result=lsim(Original_System_State_Space, U_ff, t);
plot(y_result, 'LineWidth', LineWidth);
xlabel('time(ms)');
ylabel('Voltage');


%% Close loop
K=0.1;
A_cl=(A-B*K);
B_cl=B;
C_cl=C;
D_cl=D;

input_close_loop=U_ff+K*X_ref;

system_close_loop=ss( A_cl, B_cl, C_cl, D_cl );
y_close_loop_position=lsim(system_close_loop, input_close_loop(:,1), t);

figure(5);
plot(y_close_loop_position, 'LineWidth', LineWidth);



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


