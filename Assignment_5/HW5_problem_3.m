close all; clear all;
%% defining the desired output and filtering it two times
tin = 10; tup = 11;  tdown=15; tback=16; tf = 25; delt = 0.01; ymax = 1;
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
Wf = 1; % choose a filter break frequency in Hz
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
[mm,nn] = size(y2d);
time = 0:delt:(mm-1)*delt;
U = [yd y1d];
%% Define the system G; relative degree = 1
r=1;
[A,B,C,D] = tf2ss([1 -7 2 10],[1 10 35 50 24]);
T=[C;
    0 1 0 0 ;
    0 0 1 0 ;
    0 0 0 1 ];
T_in = inv(T);
A_hat = (T*A*T_in) -((T*B*C*A*A*T_in)/(C*A*B));
B_hat = T*B/(C*A*B);  % nominal system with inverse input
A_int = A_hat(2:4,2:4); % internal dynamics 
B_int = [A_hat(2:4,1) B_hat(2:4)];

%% mdc decomposition
[As, Au, Anh, A_dec, T_mdc]=mdc(A_int, 'd');
check_A=inv(T_mdc)*A_int*T_mdc; % check_A should be eqqual to A_dec
As=A_dec(1,1);
Au=A_dec(2:3,2:3);
Bdiag_mdc=inv(T_mdc)*B_int; 
Bs=Bdiag_mdc(1,:);
Bu=Bdiag_mdc(2:3,:);

%% simulating the unstable portion  of the internal dynamics 
Utemp  = flipud(U);
% (-Bu*( [ 1;0;0] ) ) ./ (Au)
% [yu,xu] = lsim(-Au,-Bu,[1 1],[0],Utemp,time, [1.1903 0 ; 0 1.1903]  );
% -(1./Au) *Bu.*( Utemp(1,:) )
[yu,xu] = lsim(-Au,-Bu,[1 1],[0], Utemp, time, -(inv(Au)) *Bu*( Utemp(1,:)' )   );
xu = flipud(xu);
%% simulating the stable portion  of the internal dynamics 
[ys,xs] =  lsim(As,Bs,[1],[0],U,time,-(inv(As))*Bs*(U(1,:)'));
figure(2); clf; subplot(211), plot(time,xs)
xlabel('time'); ylabel('xs')
subplot(212), plot(time,xu)
xlabel('time'); ylabel('xu')
%% Calculating the input
Uff = (inv(C*B))*(y1d' -C*A*T_in*( [yd'; T_mdc*[xs';xu']])); % important
y_result_1=lsim(A,B,C,D, Uff, time);
figure(3); clf; subplot(211), plot(time,Uff)
xlabel('time'); ylabel('U_{ff}');
subplot(212), plot(time, y_result_1)
xlabel('time'); ylabel('y with u_{ff}');
%% calculate the error
for i =1:1:length(y_result_1)
    error_problem_1(i)=y_result_1(i)-yd(i);
end


%% Preview time
Tp=[1:1:5];

%% Tp-1
N=(tf-Tp(1))/delt+1;
for i=1:1:N
    new_t=[(i-1)*delt:delt:Tp(1)+(i-1)*delt];
    Utemp_new = interp1( time , Utemp, new_t);
    [yu,xu] = lsim(-Au,-Bu,[1 1],[0], Utemp_new, new_t, -(inv(Au)) *Bu*( Utemp(1,:)'));
    xu = flipud(xu);
    eta_1(i,:)=xu(1,:);
end
    % stable part
[ys,xs] =  lsim(As,Bs,[1],[0],U,time,-(inv(As))*Bs*(U(1,:)'));
xs = xs(1:length(eta_1));
y1d =y1d(1:length(eta_1));
yd=yd(1:length(eta_1));
    % inverse command
Uff = (inv(C*B))*(y1d' -C*A*T_in*( [yd'; T_mdc*[xs'; eta_1']]));
Uff=Uff(1:length(new_t));    
    % calculate the error
y_result_2_1=lsim(A,B,C,D, Uff, new_t);
for i =1:1:length(y_result_2_1)
    error_problem_2_1(i)=y_result_2_1(i)-yd(i);
end


%% Tp-2
N=(tf-Tp(2))/delt+1;
for i=1:1:N
    new_t=[(i-1)*delt:delt:Tp(2)+(i-1)*delt];
    Utemp_new = interp1( time , Utemp, new_t);
    [yu,xu] = lsim(-Au,-Bu,[1 1],[0], Utemp_new, new_t, -(inv(Au)) *Bu*( Utemp(1,:)'));
    xu = flipud(xu);
    eta_2(i,:)=xu(1,:);
end
    % stable part
[ys,xs] =  lsim(As,Bs,[1],[0],U,time,-(inv(As))*Bs*(U(1,:)'));
xs = xs(1:length(eta_2));
y1d =y1d(1:length(eta_2));
yd=yd(1:length(eta_2));
    % inverse command
Uff = (inv(C*B))*(y1d' -C*A*T_in*( [yd'; T_mdc*[xs'; eta_2']]));
Uff=Uff(1:length(new_t));    
    % calculate the error
y_result_2_2=lsim(A,B,C,D, Uff, new_t);
for i =1:1:length(y_result_2_2)
    error_problem_2_2(i)=y_result_2_2(i)-yd(i);
end

%% Tp-3
N=(tf-Tp(3))/delt+1;
for i=1:1:N
    new_t=[(i-1)*delt:delt:Tp(3)+(i-1)*delt];
    Utemp_new = interp1( time , Utemp, new_t);
    [yu,xu] = lsim(-Au,-Bu,[1 1],[0], Utemp_new, new_t, -(inv(Au)) *Bu*( Utemp(1,:)'));
    xu = flipud(xu);
    eta_3(i,:)=xu(1,:);
end
    % stable part
[ys,xs] =  lsim(As,Bs,[1],[0],U,time,-(inv(As))*Bs*(U(1,:)'));
xs = xs(1:length(eta_3));
y1d =y1d(1:length(eta_3));
yd=yd(1:length(eta_3));
    % inverse command
Uff = (inv(C*B))*(y1d' -C*A*T_in*( [yd'; T_mdc*[xs'; eta_3']]));
Uff=Uff(1:length(new_t));    
    % calculate the error
y_result_2_3=lsim(A,B,C,D, Uff, new_t);
for i =1:1:length(y_result_2_3)
    error_problem_2_3(i)=y_result_2_3(i)-yd(i);
end


%% TP-4
N=(tf-Tp(4))/delt+1;
for i=1:1:N
    new_t=[(i-1)*delt:delt:Tp(4)+(i-1)*delt];
    Utemp_new = interp1( time , Utemp, new_t);
    [yu,xu] = lsim(-Au,-Bu,[1 1],[0], Utemp_new, new_t, -(inv(Au)) *Bu*( Utemp(1,:)'));
    xu = flipud(xu);
    eta_4(i,:)=xu(1,:);
end
    % stable part
[ys,xs] =  lsim(As,Bs,[1],[0],U,time,-(inv(As))*Bs*(U(1,:)'));
xs = xs(1:length(eta_4));
y1d =y1d(1:length(eta_4));
yd=yd(1:length(eta_4));
    % inverse command
Uff = (inv(C*B))*(y1d' -C*A*T_in*( [yd'; T_mdc*[xs'; eta_4']]));
Uff=Uff(1:length(new_t));    
    % calculate the error
y_result_2_4=lsim(A,B,C,D, Uff, new_t);
for i =1:1:length(y_result_2_4)
    error_problem_2_4(i)=y_result_2_4(i)-yd(i);
end

%% Tp-5
N=(tf-Tp(5))/delt+1;
for i=1:1:N
    new_t=[(i-1)*delt:delt:Tp(5)+(i-1)*delt];
    Utemp_new = interp1( time , Utemp, new_t);
    [yu,xu] = lsim(-Au,-Bu,[1 1],[0], Utemp_new, new_t, -(inv(Au)) *Bu*( Utemp(1,:)'));
    xu = flipud(xu);
    eta_5(i,:)=xu(1,:);
end
%% stable part
[ys,xs] =  lsim(As,Bs,[1],[0],U,time,-(inv(As))*Bs*(U(1,:)'));
xs = xs(1:length(eta_5));
y1d =y1d(1:length(eta_5));
yd=yd(1:length(eta_5));
%% inverse command
Uff = (inv(C*B))*(y1d' -C*A*T_in*( [yd'; T_mdc*[ xs'; eta_5']]));
Uff=Uff(1:length(new_t));    
 %% calculate the error
y_result_2_5=lsim(A,B,C,D, Uff, new_t);
for i =1:1:length(y_result_2_5)
    error_problem_2_5(i)=y_result_2_5(i)-yd(i);
end



