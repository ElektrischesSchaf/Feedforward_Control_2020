clc;
clear;
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
t=0:0.1:10;
x1=1*[t>=0];
x2=-2*[t>=1];
x3=2*[t>=3];
x4=-1*[t>=4];
U_ff=x1+x2+x3+x4;
% plot(t,input); % stem(n,x);
% title('Unit Step Signal - Cont');

%% Inverse System
Inverse_System_State_Space=ss(A_inv, B_inv, C_inv, D_inv);
%% 
yd=lsim(Inverse_System_State_Space, U_ff, t);
plot(U_ff);
hold on
plot(yd);

%%
y_position=lsim(Original_System_State_Space, yd, t);
hold on
plot(y_position);
y_position=100.*y_position;

%% Differentiate y_position twice into y_acceleration
for i=1:length(t)-1
    y_velocity(i)=y_position(i+1)-y_position(i);
end

for i=1:length(t)-2
    y_acceleration(i)=y_velocity(i+1)-y_velocity(i);
end
% clf;
plot(y_acceleration);
legend('yd','voltage input','y position', 'y acceleration');
