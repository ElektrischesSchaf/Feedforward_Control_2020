% this code show how to use optimal inverse for tracking triangular wave
%
% Szu-Chi Tien Apr. 30, 2008
%
close all;
clear all;

%====================================================
% assume we get a system model, sys,  from experiment
sys=tf([1000000],[1 500 1000000]);
figure(1)
bode(sys)

%===================================================





%==================================================
% set desired trajectory
fs = 1000;   %1Khz sampling
dt=1/fs;
t = 0:dt:0.1-dt;
y_temp = sawtooth(2*pi*25*t, 0.5); % 25Hz triangular wave 
y_temp=y_temp(11:90);
t=0:dt:1.1-dt;
yd=[zeros(1,500), y_temp,zeros(1,500)];
t=0:dt:(length(yd)-1)*dt;
figure(2)
plot(t, yd); axis([ 0 1.1 -1 1]);
xlabel('time(s)')
ylabel('yd')
legend('25Hz triangular wave')
%=================================================



%==================================================================
%	Set optimal inversion weightings for states and frequency
scale_q = 1;    % increase to improve tracking precision 
scale_r = 1;    % increase to reduce the magnitude of inputs

omega     = [0    25   50     55   80    100   125  130    200  ]*2*pi;
q1  =       [1    1     1     0     0      0     0    0      0  ]*scale_q;
r1  =       [0    0     0     1     1      1     1    1      1  ]*scale_r;


q2  =       [1    1     1     1     1      1     1    0      0  ]*scale_q;
r2  =       [0    0     0     0     0      0     0    1      1  ]*scale_r;
%========================================================================




%========================================================================
% find the optimal input and optimal output
[Z, P, K]=zpkdata(sys,'v');
[uopt1, yopt1]=optinv_s_tf(K*poly(Z), poly(P),r1,q1,omega,yd,t);  % tracking up to 80Hz
[uopt2, yopt2]=optinv_s_tf(K*poly(Z), poly(P),r2,q2,omega,yd, t); % tracking up to 130Hz

figure(3)
subplot(2,1,1)
plot(omega/2/pi, q1, omega/2/pi, r1)
xlabel('frequency (Hz)')
ylabel('Q and R')
legend('q1','r1')
title('Q and R setting for uopt1')
axis([0 200 -0.5 1.5]);

subplot(2,1,2)
plot(omega/2/pi, q2, omega/2/pi, r2)
xlabel('frequency (Hz)')
ylabel('Q and R')
legend('q2','r2')
title('Q and R setting for uopt2')
axis([0 200 -0.5 1.5]);


figure(4)
subplot(2,1,1)
plot(t, yd, t, yopt1, t, yopt2)
legend('yd', 'yopt1', 'yopt2');

subplot(2,1,2)
plot(t, uopt1, 'g', t, uopt2,'r')
legend('uopt1', 'uopt2');
xlabel('time(s)')