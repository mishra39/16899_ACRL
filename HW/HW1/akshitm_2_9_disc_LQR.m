close all
clear all
clc

%% Parameters
syms px_k py_k vel_k theta_k dv_k dtheta_k real
lambda_k1 = sym('lam_k1_',[1 4]);
lambda_k = sym('lam_k_',[1 4]);
x0 = [0;0;0;0];% Intial State
xg = [10;0;0;0]; % Final Goal state
Q = eye(4);
Qf = 0.05;% Final Cost
R = eye(2);
S = Qf*eye(4);
N = 10;

A = [0 0 cos(x0(3)) 0;
    0 0 sin(x0(3)) 0.00001;
    0 0 0 0;
    0 0 0 0];
B = [0 0;
    0 0;
    0.01  0;
    0 0.01];


lambda = cell(N,1);
f = [vel_k*cos(theta_k)*dt;
    vel_k*sin(theta_k)*dt;
    dv_k*dt;
    dtheta_k*dt];
%%
for N = 10:1:1
    lam_n_1 = S*(x_n_1-xg);
    u_n_1 = -inv(R)*B'*lam_n_1;
    lam_k = Q*(x_n_1-xg) + 
end
l_k = 0.5*((xk-xg)'*Q*(xk-xg));
l_N = 0.5*((xk-xg)'*S*(xk-xg));
lambda_k1 = l_N;

f = [vel_k*cos(theta_k)*dt;
    vel_k*sin(theta_k)*dt;
    dv_k*dt;
    dtheta_k*dt];

f_xk_uk = f + xk;
x0_opt = x0;
xk = subs(xk,x0);
for k = 0:N
    df = [-1 0 cos(theta_k)*dt -vel_k*sin(theta_k)*dt;
           0 -1 sin(theta_k)*dt vel_k*cos(theta_k)*dt;
           0 0 1 0;
           0 0 0 1];
    lambda_k = ((xk - xg)*Q) + df'*lambda_k1;
    H_k = l_k + lambda_k1*f_xk_uk;
    
end

% u = -inv(R+B'*P_k1*B)*(B'*P_k1*A*x);