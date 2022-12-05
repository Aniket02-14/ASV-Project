clc
n=40000;
v=10;
KP=1;
KD=10;
KI=0.01;
delm=45*(pi/180);
T=107.3;
K=0.185;
x_limit=n;
y_limit=n;

chii=zeros(1,n+1);
ri=zeros(1,n+1);
xi=zeros(1,n+1);
yi=zeros(1,n+1);
%del=35*(pi/180);
h=0.1;%%%time step
%%(t*d(r)+r=k*del)
t=0:h:n*h;
for i=1:n
    if(i==1)
    KII=0;
    else
    KII=KI*trapz(t(1:i),(delm-chii(1:i)));
    end 

    f1=@(t,ri,chii) (K*KP*delm-(1+K*KD)*ri-(K*KP*chii)+K*KII)/T;


    k1y1 = h*f1(t(i),ri(i),chii(i));

    k2y1 = h*f1(t(i)+h/2,ri(i)+(k1y1/2),chii(i)+(k1y1/2));

    k3y1 = h*f1(t(i)+h/2,ri(i)+(k2y1/2),chii(i)+(k2y1/2));

    k4y1 = h*f1(t(i)+h,  ri(i)+k3y1,chii(i)+k3y1);
 

   ri(i+1)=ri(i) + (k1y1 + 2*k2y1 + 2*k3y1 +k4y1)/6;
   chii(i+1)=chii(i)+ri(i+1)*h;
    xi(i+1)=xi(i)+v*cos(chii(i+1))*h;
    yi(i+1)=yi(i)+v*sin(chii(i+1))*h;
end
syms x
y=tan(delm)*x;
figure
%tt=tiledlayout(2,1);
%nexttile
fplot(y);
hold on
plot(yi,xi)
hold off
axis([0 x_limit 0 y_limit]);
title("TRAJACTORY");
legend("y=m*x","TRAJ");



figure
tr=tiledlayout(2,2);
nexttile
plot(t,ri)
title("YAW RATE");

nexttile
yline(delm*(180/pi));
hold on
plot(t,chii.*(180/pi))
hold off
title("CHII(YAW_ANGLE) in DEG");

nexttile
plot(t,xi)
%hold on
title("CO-ORDINATES POSITIONS X-(NORTH)");
%legend("X-(NORTH)","Y-(EAST)");

nexttile
plot(t,yi)
%hold off
title("CO-ORDINATES POSITIONS Y-(EAST)");
%legend("X-(NORTH)","Y-(EAST)");





