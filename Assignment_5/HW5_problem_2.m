%% Hectors model
close all; clear all;
m=10; 
M=[10 0 0 ; 0 0 0 ; 0 0 12]+(1/12).*[2 1 0;1 4 1;0 1 2];
Bh=[1;0;0];
Ke=2.*[1 -1 0;-1 2 -1;0 -1 1];
Ce=0.1.*Ke;
A = [
    0 0 0 1 0 0 ;
    0 0 0 0 1 0 ;
    0 0 0 0 0 1 ;
    -inv(M)*[Ke Ce];
    ];
B = [0 ; 0 ; 0 ; inv(M)*Bh];
C = [0 0 1 0 0 0];
D = [0];
[num,den] = ss2tf(A,B,C,D);
r=2 ; % relative degree is 2
Sys = ss(A,B,C,D);

fprintf("tzeros(Sys)");
tzero(Sys) % finding the zeros of the system


%% internal dynamics
T = [
    C;
    C*A^(r-1);
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    ];
T_in = inv(T);
A_hat = (T*A*T_in) -((T*B*C*A*A*T_in)/(C*A*B));
B_hat = T*B/(C*A*B);  % nominal system with inverse input
A_int = A_hat(3:6,3:6); % internal dynamics 
B_int = [A_hat(3:6,1:2) B_hat(3:6)];

%% mdc decomposition
[As, Au, Anh, A_dec, T_mdc]=mdc(A_int, 'd');
check_A=inv(T_mdc)*A_int*T_mdc; % check_A should be eqqual to A_dec

As=A_dec(1:2,1:2);
Au=A_dec(3:4,3:4);

%% TODO problem
%  As = [-3.84   0.01; 
%         0.01 -3.84];
% Au = [6.24  0.01; 
%       0.01 6.24];
%%
  
Bdiag_mdc=inv(T_mdc)*B_int; 
Bs=Bdiag_mdc(1:2,:);
Bu=Bdiag_mdc(3:4,:);

%% defining the desired output and filtering it two times
tin = 5; tup = 15;  tf = 25; delt = 0.001; ymax = 1;
ramp = ymax/(tup-tin);
t1 = 0:delt:tin; 
t2=max(t1)+delt:delt:tup;
t3=max(t2)+delt:delt:tf;
y0 = zeros(size(t1));
y1 = ramp*(t2-max(t1));
y2 = max(y1)*ones(size(t3));
t = 0:delt:tf; y = [y0 y1 y2];
figure(1); clf; subplot(211), plot(t,y); 
axis([0,tf,-0.1,1.2*ymax])
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
% numerical differentiation
y1d = diff(yd,1)/(delt^1); % the first derivative of yd
y1d = [y1d;0];
y2d = diff(yd,2)/(delt^2); % the second derivative of yd
y2d = [y2d;0;0];
[mm,nn] = size(y2d);
time = 0:delt:(mm-1)*delt;
figure(1); clf; 
subplot(311), plot(time,yd);xlabel('time'); ylabel('yd')
subplot(312), plot(time,y1d);xlabel('time'); ylabel('yd^1')
subplot(313), plot(time,y2d);xlabel('time'); ylabel('yd^2')

U = [yd y1d y2d];
%% simulating the unstable portion  of the internal dynamics 
Utemp  = flipud(U);
% (-Bu*( [ 1;0;0] ) ) ./ (Au)
% [yu,xu] = lsim(-Au,-Bu,[1 1],[0],Utemp,time, [1.1903 0 ; 0 1.1903]  );
[yu,xu] = lsim(-Au,-Bu,[1 1],[0],Utemp,time, -(1./Au) *Bu*( Utemp(1,:)' )   );
xu = flipud(xu);

%% simulating the stable portion  of the internal dynamics 
[ys,xs] =  lsim(As,Bs,[1 1],[0],U,time,-(1./As)*Bs*(U(1,:)'));
figure(2); clf; subplot(211), plot(time,xs)
xlabel('time'); ylabel('xs')
subplot(212), plot(time,xu)
xlabel('time'); ylabel('xu')

%% Calculating the input
Uff = (1./(C*A*B))*(y2d' -C*A*A*T_in*( [yd';y1d'; T_mdc*[xu';xs']]));
figure(3); clf; subplot(211), plot(time,Uff)
xlabel('time'); ylabel('U_{ff}')

%% Verify with forward simulation with feedback
Kp = 0.1; Kd = 0.1;
Acl = A -B*C*Kp -Kd*B*C*A;
eig(Acl)
Ut = Uff +Kp*yd' +Kd*y1d';
[yf,xf]=lsim(Acl,B,C,D,Ut,time);
subplot(212), plot( time, yd, '*', time, yf, 'g' )
xlabel('time'); ylabel('ouput achieved in green');
legend('yd', 'y with uff+ufb')