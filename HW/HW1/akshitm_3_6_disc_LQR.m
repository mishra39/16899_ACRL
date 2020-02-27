close all
clear all
clc

syms px_k py_k vel_k theta_k dv_k dtheta_k real
%% Parameters
x0 = [0;0;0;0];% Intial State
xg = [10;0;0;0]; % Final Goal state

uk = [dv_k;dtheta_k];
N = 1000; %Horizon Size

% Reward/Penalty matrices for cost function
R = eye(2);
Q = eye(4);
S = eye(4);
N = 1000;
T = 10; % Final time. But this is infinite time horizon??
dt = 0.01;

%Vector size
vec_len = (T/dt) +1;

% Linearized System Dynamics
A = [0 0 1 0;
    0 0 0 0.00001;
    0 0 0 0;
    0 0 0 0];
B = [0 0;
    0 0;
    1 0;
    0 1;];

C = zeros(1,4);
D = 0;

Ts = 0.01;
sysc = ss(A,B,C,D);
sysd = c2d(sysc,Ts);
Ad = sysd.A;
Bd = sysd.B;
Bd = [0 0; 0 0;0.01 0;0 0.01];
Pk = cell(N,1);
K = cell(1,N);
xk = cell(N,1);
u = cell(N,1);
Pk{N} = S;

% Backward Pass
for k = N-1:-1:1
    Pk{k} = Q + Ad' * Pk{k+1} * Ad - Ad' * Pk{k+1} * Bd * inv(R+Bd' * Pk{k+1} * Bd) * Bd' * Pk{k+1};
    K{k} = -inv(R + Bd' * Pk{k+1} * Bd) * Bd' * Pk{k+1} * Ad;
end

xk{1} = x0;

%Forward Pass
for k = 1:1:N-1
    u{k} = K{k}*(xk{k}-xg);
    xk{k+1} = Ad*xk{k} + Bd*u{k};
end
px = []
py = []
thta = []
v  = []
for k = 1:1:N
    px(k) = (xk{k}(1));
    py(k) = (xk{k}(2));
    v(k) = (xk{k}(3));
    theta(k) = (xk{k}(4));
end
t = 0:.01:9.99;
% State Plots
figure(2)
plot(t,px,'-o',t,py,'-x', t,v,'-x', t,theta,'-o')
title('State Plots for finite time LQR');
xlabel('Time t');
ylabel('State Variables');
legend('Position in x','Position in y', 'Velocity', 'theta','Location','SouthEast')

% Co-state plot

for k = 1:1:N
    lambda_k(k:k+3) = Pk{k}*(xk{k} - xg);
end
figure(2)
plot(t,lambda_k(1:1000),'-x')
title('Co-State Plot for finite time LQR');
xlabel('Time t');
ylabel('Co-State Variables');

%% Control Trajectory plot
vel_kk =[];
theta_kk = [];
for k = 1:1:N-1
    vel_kk(k)= u{k}(1);
    theta_kk(k) = u{k}(2);
end


figure(3)
plot(t(1:999),vel_kk,'-o',t(1:999),theta_kk,'-x')
title('Control Input Plot for finite time LQR');
xlabel('Time t');
ylabel('Control Inputs');
legend('Velocity Derivative Input','Theta dot input','SouthEast')
%% 3.5
% l_k = 0.5*((xk-xg)'*Q*(xk-xg));
% l_N = 0.5*((xk-xg)'*S*(xk-xg));
% lambda_k1 = l_N;
% f = [vel_k*cos(theta_k)*dt;
%     vel_k*sin(theta_k)*dt;
%     dv_k*dt;
%     dtheta_k*dt];
% 
% f_xk_uk = f + xk;
% x0_opt = x0;
% xk = subs(xk,x0);
% for k = 0:N
%     df = [-1 0 cos(theta_k)*dt -vel_k*sin(theta_k)*dt;
%            0 -1 sin(theta_k)*dt vel_k*cos(theta_k)*dt;
%            0 0 1 0;
%            0 0 0 1];
%     lambda_k = ((xk - xg)*Q) + df'*lambda_k1;
%     H_k = l_k + lambda_k1*f_xk_uk;
%     
% end
% 
% % u = -inv(R+B'*P_k1*B)*(B'*P_k1*A*x);