clc; clear;
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6]);
time=[0:0.1:8];
uff=[1:0.1:9];
U=[time; uff];

t_start=0; t_end=8;
tspan=[t_start t_end];
IC=[0 0 0];
[t, x_n]=ode45(@(t, x_n) fun(t, x_n, U), tspan, IC, options);
% note: time and t are different

% function xdiff= fun(t, x_n, flag, U)
function xdiff= fun(t, x_n, U)
    J=1;

    time=U(1,:);
    uff=U(2,:);
    u = interp1(time, uff, t);
%     uff = @(time_handle) 50*cos(time_handle); % mine
    
    xdiff=zeros(3,1);
    xdiff(1)=x_n(2);
    xdiff(2)=(x_n(3)-sin(x_n(1))-x_n(2))/J;
    
    xdiff(3)=-x_n(3)-x_n(2)+u;
%     xdiff(3)=-x_n(3)-x_n(2)+uff(t); % mine

end