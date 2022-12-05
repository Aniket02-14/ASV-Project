clc
clear all
close all
T = 107.3;
k = 0.185;
kp = 1;
kd = 10;
ki = 0.01;
v = 10;
del_t = 0.1;
psi_d = 45*pi/180;
%t = 0:del_t:4000;
a = 1+k*kd;
b = k*kp;
c = k*ki;
d = psi_d*k*ki;
%% Example 4
% Solve y'''(t)+4y''(t)+6y'(t)+4y(t)=gamma(t),  y''(0)=0, y'(0)=-1, y(0)=0
q0 = [0; 0; 0];         % Initial Condition (vector)
h=0.1;                   % Time step
t = 0:h:5000;               % t goes from 0 to 2 seconds.
A = [0 1 0; 0 0 1; -a -b -c];  % A Matrix
B = [0; 0; 1];                 % B Matrix


qstar = zeros(3,length(t));  % Preallocate array (good coding practice)

qstar(:,1) = q0;             % Initial condition gives solution at t=0.
for i=1:(length(t)-1)
  k1 = A*qstar(:,i)+B*d;            % Approx for y gives approx for deriv
  q1 = qstar(:,i)+k1*h/2;           % Intermediate value
  k2 = A*q1+B*d;                    % Approx deriv at intermediate value.
  qstar(:,i+1) = qstar(:,i) + k2*h; % Approx solution at next value of q
end

i = 1;
x = [0];
while i < length(t)
    x(end+1) = x(end) + v*sin(qstar(1,i))*del_t;
    i = i+1;
end
y = [0];
for k = 1:length(x)
    if k == 1
        y(k) = 0 + v*cos(qstar(1,k))*del_t;
    else
        y(k) = y(k-1) + v*cos(qstar(1,k))*del_t;
    end
end
plot(x,y,x,x)
%xlabel("x")
%ylabel("y")
%title("x vs. y")