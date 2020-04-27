% =========================================================
% this function is saved as ¡¥findx_n¡¦ in the same folder of your main program 
% ========================================================

function xdiff= fun(t, x_n, flag, J, U)
time = U(:,1); 
uff = U(:,2); 
u = interp1(time, uff, t); 
xdiff=zeros(3,1);
xdiff(1)=x_n(2);
xdiff(2)=(x_n(3)-sin(x_n(1))-x_n(2))/J;
xdiff(3)=-x_n(3)-x_n(2)+u;
return
