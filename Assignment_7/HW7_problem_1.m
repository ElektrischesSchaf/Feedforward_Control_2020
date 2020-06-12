clc; clear;
% num=[1.0304 0.3770];
% den=[1.0 -0.2233 0.0561];
% h=dimpulse(num,den,10);
% figure(1);
% stem(0:1/2000:9*1/2000,h);
% xlabel('time(sec)');
% ylabel('impulse response');

G_hat_s_domain=tf([1],[0.8 8 30]);
Hd=c2d(G_hat_s_domain,0.2, 'zoh')
step(G_hat_s_domain, '-', Hd,'--');


num=[ 0.01274 0.006455];
den=[ 2 -0.5594 0.1353];
h=dimpulse(num,den,10);
figure(1);
stem( 0:1/5:9*1/5, h );
xlabel('time(sec)');
ylabel('impulse response');

t = 0:0.05:10;
yd = 0.1*sin(2*pi*0.75*t);
yn = 0.005*rand(1, length(t) );
