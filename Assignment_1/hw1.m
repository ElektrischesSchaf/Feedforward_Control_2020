clear;
%% initialize the system
numerator=[11.88 4.977 539.6 129.9 5625]
denominator=[1 1.169 50.3 45.94 685.3 391.7 1952]
sys=tf(numerator, denominator);

%% Bode Plot of the system
% bode(sys); % change frequency unit to KHz in "View" -> "Property Editor"

%% Step response of the system
% plot(step(sys));

%% The relative degree of the system
% relative degree =  6-4 =2

%% Original System State Space
csys= canon(sys,'companion')
A=csys.A;
B=csys.B;
C=csys.C;
D=csys.D;

%% inverse system State Space
A_inv=A-((B*C*A*A)/(C*A*B));
B_inv=(B)/(C*A*B);
C_inv=(-C*A*A)/(C*A*B);
D_inv=1/(C*A*B);

%% signal generation
% https://electrosome.com/signal-generation-in-matlab/
n1=0; n2=10;
t=n1:0.1:n2;
x1=1*[t>=0];
x2=-2*[t>=1];
x3=2*[t>=3];
x4=-1*[t>=4];
input=x1+x2+x3+x4;
% plot(t,input); % stem(n,x);
% title('Unit Step Signal - Cont');

inverse_system=ss(A_inv, B_inv, C_inv, D_inv);
lsim(inverse_system, input, t);