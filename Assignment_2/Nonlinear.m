clc; clear;
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6]);

% U=[time; uff];

[time, state]=ode45(@fun,[0 8],[0 0 0],options);
% note: time and t are different

function xdiff= fun(t, x_n, flag)
    J=1;
    uff = @(time) 50*cos(time);
%   u = interp1(time, uff, t);
    xdiff=zeros(3,1);
    xdiff(1)=x_n(2);
    xdiff(2)=(x_n(3)-sin(x_n(1))-x_n(2))/J;
    xdiff(3)=-x_n(3)-x_n(2)+uff(t);
end