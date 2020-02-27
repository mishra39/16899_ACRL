close all
clear all
clc
%% Parameters
x0 = [0 ;0; 0; 0];% Intial State
xg = [1 ;1 ;0.00001 ;0]; % Final Goal state
N = 10; %Horizon Size

% Reward/Penalty matrices for cost function
R = 1;
Q = [1 0 0 0; 0 1 0 0; 0 0 1 0;0 0 0 1];

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
C0 = rank(ctrb(A,B));
disp('The rank of controllability matrix is')
disp(C0)
disp('The linearized system is uncontrollable')
tspan = 0:0.01:10; % time span for ode solver

K = lqr(A,B,Q,R)

% Solve ARE to find update state x
[t x] = ode45(@(t,x)func(t, x,-K*(x-xg)), tspan, x0) ;

plot(t,x(:,1),'-o',t,x(:,2),'-x', t,x(:,3),'-x', t,x(:,4),'-o')
title('Solution of inifinite time LQR');
xlabel('Time t');
ylabel('State Variables');
legend('Position in x','Position in y', 'Velocity', 'theta')

% Apply Optimal control lab

% P = icare(A,B,Q,R,[],[],G)
% P_dot_t = P_t*A_t + A_t'.P_t - P_t*B_t*R_t
% Solve for S backwards in time using Ricatti
% [t_out_S, S_ode] = ode45(@(t,S) Ricatti(S), [
% 
% 
% function dXdt = Riccati(t,X,A,B,Q)
% 
% X = reshape(X,size(A));
% dXdt = A.'*X + X*A - X*B*B.'*X + Q;
% dXdt = -dXdt(:);
% end
% A = [0 0 cos(theta_r) -v_r*sin(theta_r);
%     0 0 sin(theta_r) v_r*cos(theta_r);
%     0 0 0 0;
%     0 0 0 0];
