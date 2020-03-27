clc;
clear;
L=1; M=1; m=1; g=1; B=1;
syms eta1 eta2;
syms Xc Xc_d theda theda_d;
z=[Xc Xc_d theda theda_d];
% [Xc_dd theda_dd]=[M+m m*L*cos(theda) ; m*L*cos(theda) m*L^2]'*[m*L*sin(theda)*]
