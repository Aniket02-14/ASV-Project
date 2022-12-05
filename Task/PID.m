clc
clear
close all
T = 107.3;
k = 0.185;
kp = 1;
kd = 10;
ki = 0.01;
v = 10;
del_t = 0.1;
psi_d = 45*pi/180;
t = 0:del_t:300;
a = (1+k*kd)/T;
b = k*kp/T;
c = k*ki/T;
e = psi_d;
i = e*del_t;
d = psi_d*k*kp;
psi = zeros(length(t), 1);
psi(1,1) = 0;
r = zeros(length(t),1);
r(1,1) = 0;
x = zeros(length(t),1);
y = zeros(length(t),1);
x(1,1) = 0;
y(1,1) = 0;
for j = 2 :length(t)
    k1 = d-b*e+c*i-a*r(j-1,1);
    k2 = d-b*(e-r(j-1,1)*del_t/2)+c*(i-e*del_t/2)-a*(r(j-1,1)+del_t*k1/2);
    k3 = d-b*(e-r(j-1,1)*del_t/2)+c*(i-e*del_t/2)-a*(r(j-1,1)+del_t*k2/2);
    k4 = d-b*(e-r(j-1,1)*del_t)+c*(i-e*del_t)-a*(r(j-1,1)+del_t*k3);
    k_o = (k1+2*k2+2*k3+k4)/6;
    r(j,1) = r(j-1,1) + k_o;
    e = e - r(j,1)*del_t;
    psi(j,1) = psi(j-1,1) + r(j,1)*del_t;
    i = i + e*del_t;
    x(j,1) = x(j-1,1) + v*sin(psi(j,1))*del_t;
    y(j,1) = y(j-1,1) + v*cos(psi(j,1))*del_t;
    j = j+1;
end
plot(x,y,x,x)
