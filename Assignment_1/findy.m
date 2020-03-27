function dy=findy(time,y);
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

