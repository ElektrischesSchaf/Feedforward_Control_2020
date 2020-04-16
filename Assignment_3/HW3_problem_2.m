clc; clear;
r=2;
M=1; m=1; L=1; B=1; g=1;
%% defining the desired output (and filtering it twice)
tin = 1; tup =3;  tf = 8; delt = 0.001; ymax = 10;

ramp = ymax/(tup-tin);
t1 = 0:delt:tin; 
t2=max(t1)+delt:delt:tup;
t3=max(t2)+delt:delt:tf;
y0 = zeros(size(t1));
y1 = ramp*(t2-max(t1));
y2 = max(y1)*ones(size(t3));
t = 0:delt:tf; y = [y0 y1 y2];

%% filtering the desired signal 
Wf = 1; % choose a filter break frequency in Hz
num = [Wf*2*pi]; den = [1 (Wf*2*pi)]; % first order filter
[Af,Bf,Cf,Df] = tf2ss(num,den);
Sys_f = ss(Af,Bf,Cf,Df); % first order filter system
Sysf = Sys_f*Sys_f*Sys_f*Sys_f*Sys_f; % fifth order filter
[yd,xtemp]= lsim(Sysf,y,t); % filter the desired signal
figure(1);
subplot(211), plot(t,y); 
axis([0 8 -1 11]);
xlabel('time(s)'); ylabel('y (unfiltered)')
subplot(212), plot(t,yd); 
xlabel('time(s)'); ylabel('yd (filtered)')
axis([0 8 -1 11]);

%% y_desired_double_dot
for i=1:length(y)-1
    y_velocity(i)=yd(i+1)-yd(i);
end
cc=[0]
y_velocity=[y_velocity cc];
for i=1:length(y_velocity)-1
    y_acc(i)=y_velocity(i+1)-y_velocity(i);
end
y_acc=[y_acc cc];

%% Calculate inveres response
% https://stackoverflow.com/questions/26916066/variable-input-changing-in-time-ode45-matlab
eta1=0; eta2=0;
IC=[eta1, eta2];
% time_span=[0 tf]; % If with this, the response will be downsampled.
time_span=t;
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]); 
options2 = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6]); 
[time, eta] = ode45( @(time, eta) func(time, eta, y_acc, M, m, L, B, g), time_span, IC, options);

M=1; m=1; L=1; B=1; g=9.8;
ttt= [0: delt : 8];
y_desired_dot_dot=interp1( ttt, y_acc, time );
U_inv=-m*L*sin(eta(:,1)).*eta(:,2).^2+(M+m).*y_desired_dot_dot+m*L.*cos(eta(:,1)).*((-1/m*L^2)*(m*g*L.*sin(eta(:,1))+B.*eta(:,2)+m*L.*cos(eta(:,1)).*y_desired_dot_dot ));


% yee=0:0.2:8;
% U_inv=interp1( 1:41, U_inv, 0:delt:8);
figure(2);
plot(U_inv);
title('U_f_f');
U=[U_inv ,  t'];

figure(3);
[t, x_n]=ode45(@(t, x_n) cart(t, x_n, U, M, m, L, M, g), time_span, [0 0 0 0], options2);
plot(x_n(:,1));
title('Origianl system output after applied U_f_f');

%% close loop control 
% TODO find A, B, C, D
% https://www.mathworks.com/help/control/getstart/pole-placement.html
% https://www.mathworks.com/help/control/ref/place.html
new_A=[0,0,0,0;
    0 0 0 m*L*B;
    0 0 0 0;
    0 0 0 (M+m)*-B]./(m*M*L);

new_B=[0 ; m*L*L ; 0 ; -m*L];
new_C=[0 0 1 0];
new_D=0;
p=[-1, -2, -3, -4];
K=1;
% K=place(new_A, new_B, p);
A_cl=(new_A-new_B*K);
B_cl=new_B;
C_cl=new_C;
D_cl=new_D;
input_close_loop=U_inv + K*X_ref;
system_close_loop=ss(A_cl, B_cl, C_cl, D_cl);
input_close_loop=input_close_loop';
t=t';
X_ref_transpose=X_ref';
y_close_loop=lsim(system_close_loop, input_close_loop(:,1), t); % input_close_loop(:,1): the first column of close loop input



function eta_dot=func(time, eta, y_acc, M, m, L, B, g)
    
    ttt= [0: 0.001 : 8];
    y_desired_dot_dot=interp1( ttt, y_acc, time );
    
    eta_dot=zeros(2,1);
    
    eta_dot(1)=eta(2);
    eta_dot(2)=( -1/(m*L*L) )*(m*g*sin(eta(1))+B*eta(2)) + ((-1/m*L*L))*(m*L*cos(eta(1)))*y_desired_dot_dot;
end


function x_diff = cart(t, x_n , U, M, m, L, B, g)
    
    F=U(:,1);
    time=U(:,2);    
    u = interp1(time, F, t);
    
    x_diff=zeros(2,1);
    x_diff(1)=x_n(2);
    x_diff(3)=x_n(4);
    A_s=[m*L*L  , -m*L.*cos(x_n(3));
        -m*L.*cos(x_n(3))   ,  M+m]
    A_s=A_s.*(1/(m*M*L*L+m*m*L*L.*sin(x_n(3)).*sin(x_n(3)) ) )
    B_s=[m*L.*sin(x_n(4)).*x_n(4)^2+u;
        -m*g*L.*sin(x_n(3))-B.*x_n(4)];
    C_s=A_s*B_s;
    x_diff(2)=C_s(1);
    x_diff(4)=C_s(2);

end

