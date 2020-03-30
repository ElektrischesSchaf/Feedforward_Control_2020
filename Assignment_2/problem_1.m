clc; clear;

L=1; R=1; K=1; m=1; g=1; J=1;B=1;
A=[
    0, 1, 0;
    0, -B/J, K/J;
    0, -K/L, -R/L;
    ];
B=[ 0 ; 0 ; 1/L];
C=[1, 0, 0];
D=[0];
[numerator denominator]=ss2tf(A,B,C,D);
origianl_system_tf=tf(numerator, denominator);