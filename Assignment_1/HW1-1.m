clc;
clear;
close all
r=6-4;
A=[-1.169 -50.3 -45.94 -685.3 -391.7 -1952; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;0 0 0 0 1 0];
B=[1;0;0;0;0;0];
C=[0 11.88 4.977 539.6 129.9 5625];
D=0;
%% 1-(c)
Ainv=A-((B*C*A^(r))/(C*A^(r-1)*B));
Binv=(B)/(C*A^(r-1)*B);
Cinv=-(C*A*A)/(C*A^(r-1)*B);
Dinv=1/(C*A^(r-1)*B);
%% 1-(d)
t=0:0.01:10;
y=[1*ones(1,100) -1*ones(1,200) 1*ones(1,100) 0*ones(1,601)];
[uff,Xinv]=lsim(Ainv,Binv,Cinv,Dinv,y,t);
figure(1)
subplot(2,1,1)
plot(t,y);
xlabel('time(ms)');
ylabel('y^(^2^)');
title('Desired acceleration profile');
subplot(2,1,2)
plot(t,uff);
xlabel('time(ms)');
ylabel('u_i_n_v');
title('Feedeorward input u_i_n_v');
%% 1-(e)
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]); 
[tt,yy]=ode45(@findy,[0 10],[0 0],options);
ya=interp1(tt,yy(:,2),t);
yd=interp1(tt,yy(:,1),t);
[ysimu,xref]=lsim(A,B,C,D,uff,t);
figure(2)
hold on
plot(t,ya, '*')
plot(t,yd,'-')
plot(t,ysimu,'x')





%% sub
function dy=findy(time,y)
    dy=zeros(2,1);
    dy(1)=y(2);
    if time<1
       dy(2)=1;
    elseif (time>=1 && time<3)
       dy(2)= -1;
    elseif (time>=3 && time<4)
        dy(2)=1;
    else
        dy(2)=0;
    end
end