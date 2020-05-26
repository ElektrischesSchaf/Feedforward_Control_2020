% Inversion for the VTOL aircraft
% see paper : A different look at Output Tracking: Control of a VTOL Aircraft
%             By P. Martin, S. Devasia, and B. Paden
%             Automatica, Vol. 32, No. 1, pp.101-107, 1996
%             Please email questions to devasia@u.washington.edu

% Part 1: 
% defining the desired output and filtering it two times
tin = 1; tup =3;  tf = 8; delt = 0.001; ymax = 10;
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
Wf = 1; % choose a filter break frequency in Hz
num = [Wf*2*pi]; den = [1 (Wf*2*pi)]; % first order filter
[Af,Bf,Cf,Df] = tf2ss(num,den);
Sys_f = ss(Af,Bf,Cf,Df); % first order filter system
Sysf = Sys_f*Sys_f*Sys_f*Sys_f*Sys_f; % fifth order filter
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
% in this example 
% the same desired trajectory is used in the x and y directions
Xs= 1; Ys = 1;
xd = yd*Xs; x1d = y1d*Xs; x2d = y2d*Xs;  
yd = yd*Ys; y1d = y1d*Ys; y2d = y2d*Ys;  

figure(2); clf; subplot(311), plot(time,xd)
xlabel('time'); ylabel('xd')
subplot(312), plot(time,x1d)
xlabel('time'); ylabel('xd^1')
subplot(313), plot(time,x2d)
xlabel('time'); ylabel('xd^2')

%return

% Part 2: %%%%%% Description of the System 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 1; 
epsilon = 0.1; 
g = 9.81; 
gamma =  sqrt(lambda*g/epsilon);

% Part 3: Iterative solution of the internal dynamics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initiate the iterative procedure 
eta = zeros(2,length(time))';
theta = eta(:,1);
TransM = (1/(2*gamma))*[gamma -1; gamma 1];
Aold = [0 1; gamma*gamma 0]; 
TransM*Aold*inv(TransM);  % check internal dynamics decoupling

figure(3); clf
theta_store = [];
nsteps = 30;
%% iterative procedure 
for jj = 1:1:nsteps
    p = (lambda/epsilon)*(((x2d).*cos(theta))+((y2d).*sin(theta)) +(g*sin(theta))  -g*theta);
    [zs,zstemp] = lsim(-gamma, (-1/(2*gamma)), 1, 0, p, t);
    [zu,zutemp] = lsim(-gamma, (-1/(2*gamma)), 1, 0, flipud(p), t);
    zu = flipud(zu);
    %%%% compute the new internal dynamics
    eta = (inv(TransM)*[zs zu]')';
    error = norm(theta -eta(:,1))
    theta = eta(:,1);
    theta_1d = eta(:,2);
    theta_store = [theta_store theta];
    %%%% plot theta and internal states every iteration 
    subplot(311), plot(time,zs); hold on
    xlabel('time'); ylabel('z_s')
    subplot(312), plot(time,zu)
    xlabel('time'); ylabel('z_u'); hold on
    subplot(313), plot(time,theta)
    xlabel('time'); ylabel('\theta')
    drawnow
    zoom on; hold on
end
figure(3); hold off
%%%% plot theta changes over the iterations
figure(4); clf
plot(time,theta_store,time,theta,'r')
xlabel('time'); ylabel('iterations of \theta')


% Part 4: Computing the inverse 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1 = (-(x2d).*sin(theta)) +((y2d+g).*cos(theta));
u2 = (((x2d).*cos(theta)) +((y2d+g).*sin(theta)))/(epsilon);
figure(5); clf
subplot(211), plot(time,u1); 
xlabel('time'); ylabel('u_1')
subplot(212), plot(time,u2)
xlabel('time'); ylabel('u_2');



% Part 5: Apply the inverse to the nomial system 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = [time' u1 u2 theta theta_1d]; % an array with the time, input and theta vectors
param = [lambda epsilon g  gamma] % a vector with system parameters
[t,xa]=ode45('vtoldiff',[0 10],[0;0;0;0;0;0],[],param,U);


%%% plotting the simulation results


figure(6); clf;
subplot(311), plot(t,xa(:,1),'g',time,xd,'r');
xlabel('time'); ylabel('xd in red, xactual in green')
subplot(312), plot(t,xa(:,3),'g',time,yd,'r');
xlabel('time'); ylabel('yd in red, yactual in green')
subplot(313), plot(t,xa(:,5),'g',time,theta,'r');
xlabel('time'); ylabel('desired \theta in red, \theta in green')

return


