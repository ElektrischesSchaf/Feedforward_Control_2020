clc; clear;

%%
% r=3;
% m=1; k=1; c=0.1;
% A=[0, 1, 0, 0;
%     -k/m, -c/m, k/m, c/m;
%     0, 0, 0, 1;
%     k/m, c/m, -k/m, -c/m];
% B=[0 ; 1/m ; 0; 0];
% C=[0 0 1 0];
% D=0;
% 
% Tt=[C; C*A; C*A*A];
% Tb=[0 1 0 0];
% T=[Tt; Tb];
% T_inv=inv(T);
% T_invL=T_inv(:,[1,3]);
% T_invR=T_inv(:,4);
% 
% A_inv_red=[ Tb*A*T_invL - ( (Tb*B*C*A^r*T_invR)/ (C*A^(r-1)*B) )];
% B_inv_red=[(Tb*B) / (C*A^(r-1)*B) ];



%% 
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6]);
[t, eta]=ode45(@(t, eta) func(t,eta) , [0 10], [0 0], options);

function eta_dot= func( t, eta, flag)

    r=3;
    m=1; k=1; c=0.1;
    A=[0, 1, 0, 0;
        -k/m, -c/m, k/m, c/m;
        0, 0, 0, 1;
        k/m, c/m, -k/m, -c/m];
    B=[0 ; 1/m ; 0; 0];
    C=[0 0 1 0];
    D=0;

    Tt=[C; C*A; C*A*A];
    Tb=[0 1 0 0];
    T=[Tt; Tb];
    T_inv=inv(T);
    T_invL=T_inv(:,[1,3]);
    T_invR=T_inv(:,4);

    A_inv_red=[ Tb*A*T_invL - ( (Tb*B*C*A^r*T_invR)/ (C*A^(r-1)*B) )];
    B_inv_red=[(Tb*B) / (C*A^(r-1)*B) ];

    Y_d=zeros(2,1);
    
    eta_dot=A_inv_red*eta+B_inv_red*Y_d;
end

