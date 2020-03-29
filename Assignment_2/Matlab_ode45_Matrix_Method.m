function sdot=practice(t,s)
% IC
x0=0; v0=0; theta0=0; omega0=0;
IC=[x0, theta0, v0, omega0];

% Time Span
t0=0; tf=30; tspan=[t0, tf];

[time, state_values]=ode45(@g, tspan, IC);
x=state_values(:,1);
theta=state_values(:,2);
v=state_values(:,3);
omega=state_values(:,4);

figure(1), clf
plot(time,x), xlabel('time (s)'), ylabel('displacement (m)');
title('\theta displacement vs Time');

figure(2), clf
plot(time,theta), xlabel('time (s)'), ylabel('angular displacement (rad)');
title('\theta angular displacement vs Time');
end

function  sdot=g(t,s)
g=9.8;
c=5;
k=100;
L=7;
m1=500;
m2=200;

f=@(time) 50*cos(time);

M=[(m1+m2), m2*L ; m2*L, m2*L*L];

C=[c, 0 ; 0, 0];
K=[k, 0 ; 0, m2*g*L];

ZERO=zeros(2,2);
I=eye(2,2);

A=[ZERO, I; -inv(M)*K, -inv(M)*C];
B=[ZERO; inv(M)];
F=[f(t);0];
Svec=[s(1) ;  s(2) ; s(3) ; s(4)];

sdot=A*Svec+B*F;

end


