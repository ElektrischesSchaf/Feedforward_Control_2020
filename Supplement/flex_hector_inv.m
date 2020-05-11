% Hectors model for class
close all; clear all;
M = [12 1; 1 12]; 
C = 0.1*[1 -1; -1 1];
K = 10*C;
Bh = [1;0];
A = [0 0 1 0; 0 0 0 1; -inv(M)*[K C]]; 
B = [0;0;inv(M)*Bh];
C = [0 1 0 0 ]; D = [0];
[num,den] = ss2tf(A,B,C,D);
Sys = ss(A,B,C,D);
tzero(Sys) % finding the zeros of the system

% internal dynamics
T = [C;C*A;1 0 0 0 ; 0 0 1 0];
Tin = inv(T);
Ahat = (T*A*Tin) -((T*B*C*A*A*Tin)/(C*A*B));
Bhat = T*B/(C*A*B);  %%%%% nominal system with inverse input
Aint = Ahat(3:4,3:4); %%% internal dynamics 
Bint = [Ahat(3:4,1:2) Bhat(3:4)];
eig(Aint);

% the decoupled internal dynamics
N =  [1.0000    1.0512; 1.0000   -0.9512];
Nin = inv(N);
Adiag = N *Aint*Nin;
Bdiag = N*Bint;
As = Adiag(2,2); Bs = Bdiag(2,:); %   stable component of internal dynamics 
Au = Adiag(1,1); Bu = Bdiag(1,:); %   unstable component of internal dynamics 

%return
% defining the desired output and filtering it two times
tin = 10; tup = 15;  tf = 25; delt = 0.001; ymax = 10;
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
figure(1); clf; subplot(311), plot(time,yd)
xlabel('time'); ylabel('yd')
subplot(312), plot(time,y1d)
xlabel('time'); ylabel('yd^1')
subplot(313), plot(time,y2d)
xlabel('time'); ylabel('yd^2')

%return
U = [yd y1d y2d];
% simulating the unstable portion  of the internal dynamics 
Utemp  = flipud(U);
[yu,xu]=lsim(-Au,-Bu,[1],[0],Utemp,time,-(1/Au)*Bu*(Utemp(1,:)'));
xu = flipud(xu);
% simulating the stable portion  of the internal dynamics 
[ys,xs]=  lsim(As,Bs,[1],[0],U,time,-(1/As)*Bs*(U(1,:)'));
figure(2); clf; subplot(211), plot(time,xs)
xlabel('time'); ylabel('xs')
subplot(212), plot(time,xu)
xlabel('time'); ylabel('xu')

% Calculating the input
Uff = (1/(C*A*B))*(y2d' -C*A*A*Tin*( [yd';y1d'; Nin*[xu';xs']]));
figure(3); clf; subplot(211), plot(time,Uff)
xlabel('time'); ylabel('U_{ff}')

% Verify with forward simulation with feedback
Kp = 0.1; Kd = 0.1;
Acl = A -B*C*Kp -Kd*B*C*A;
eig(Acl)
Ut = Uff +Kp*yd' +Kd*y1d';
[yf,xf]=lsim(Acl,B,C,D,Ut,time);
subplot(212), plot(time,yd,'r',time,yf,'g')
xlabel('time'); ylabel('ouput achieved in green');
legend('yd', 'y with uff+ufb')

return

% K = [-3 -2];
% AA = A - B*C*A*A/(C*A*B) +B*K*[C;C*A]/(C*A*B);
% BB = [-B*K/(C*A*B)  B/(C*A*B)];
% CC = [(-C*A*A/(C*A*B))+(K*[C;C*A])];
% DD = [-K 1];
% NewSys = ss(AA,BB,C,D);
% NewSys2 = ss(AA,BB,CC,DD);
% 
% [yy,xx]=lsim(NewSys,U,time);
% [uff,xx]=lsim(NewSys2,U,time);
% 
% figure(2); clf;
% subplot(211), plot(time,yy,'g',time,yd,'r');
% xlabel('time'); ylabel('yd in red, yactual in green')
% subplot(212), plot(time,uff);
% xlabel('time'); ylabel('total input')
% zoom on
% 
% return
% 

