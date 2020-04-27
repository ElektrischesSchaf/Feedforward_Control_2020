time=[0:0.01:8];
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6]);
uff=ones(length(time),1);
U=[time', uff];  %uff means the input
[t, state]=ode45('findx_n',[0 8],[0 0 0],options,1,U);
% note: time and t are different
