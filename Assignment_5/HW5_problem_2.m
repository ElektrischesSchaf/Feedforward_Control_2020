%% Hectors model
close all; clear all;
m=10; 
M=[10 0 0 ; 0 0 0 ; 0 0 12]+(1/12).*[2 1 0;1 4 1;0 1 2];
Bh=[1;0;0];
Ke=2.*[1 -1 0;-1 2 -1;0 -1 1];
Ce=0.1.*Ke;
A = [
    0 0 0 1 0 0 ;
    0 0 0 0 1 0 ;
    0 0 0 0 0 1 ;
    -inv(M)*[Ke Ce];
    ];
B = [0 ; 0 ; 0 ; inv(M)*Bh];
C = [0 0 1 0 0 0];
D = [0];
[num,den] = ss2tf(A,B,C,D);
r=2 ; % relative degree is 2
Sys = ss(A,B,C,D);
tzero(Sys) % finding the zeros of the system

%% internal dynamics
T = [
    C;
    C*A^(r-1);
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    ];

