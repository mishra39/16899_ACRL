close all
clear all
clc

syms px_k py_k v_k theta_k dv_k dtheta_k real
%% Parameters
x0 = [0;0;0;0];% Intial State
xg = [10;0;0;0]; % Final Goal state


% Reward/Penalty matrices for cost function
R = eye(2);
Q = eye(4);
S = eye(4);
dt = 0.01;
N = 10;

lambda = cell(N,1);
u = cell(N,1);
xk = [px_k ,py_k ,v_k ,theta_k ]';
uk = [dv_k, dtheta_k]';
f = [xk(3)*cos(xk(4))*dt;
    xk(3)*sin(xk(4))*dt;
    dv_k*dt;
    dtheta_k*dt];

f_xk_uk = f + xk;
df_xk_uk = jacobian(f_xk_uk,xk)
B = [ 0 0; 0 0 ; dt 0; 0 dt];
xk_1 = cell(N,1);
lambda = cell(N,1);
xk_opt = cell(N,1);
u = cell(N,1);
uk_opt = cell(N,1);
xk_1{N-1} = xk;

%Backward Pass
for k = N-1:-1:1
   
   lambda{k+1} =  S*(f_xk_uk-xg);% lambda_N equation
   lambda{k} = Q*(xk - xg) + df_xk_uk' * lambda{k+1};% lambda_N-1 equation
   u{k} =  (-inv(R)*B')*lambda{k+1}; % u_N-1 equation
   xk_1{k+1} = f_xk_uk;         % x_k+1 dynamics of the system
   eqn1 = lambda{k+1} - S*(f_xk_uk-xg) == 0;
   eqn2 =  lambda{k}  - Q*(xk - xg) + df_xk_uk' * lambda{k+1} == 0;
   eqn3 = u{k} -  (-inv(R)*B')*lambda{k+1} == 0;
   eqn4 = xk_1{k+1} - f_xk_uk == 0;
   sol = solve([eqn1,eqn2,eqn3,eqn4],[px_k, py_k, v_k,theta_k,dv_k,dtheta_k]);
   px_k_sol = sol.px_k;
   py_k_sol = sol.py_k;
   v_k_sol = sol.v_k;
   theta_k_sol = sol.theta_k;
   dv_k_sol = sol.dv_k;
   dtheta_k_sol = sol.dtheta_k;
   xk = [px_k_sol, py_k_sol, v_k_sol, theta_k_sol]'; 
   
   f = [xk_1{k+1}(3)*cos(xk_1{k+1}(4))*dt;
        xk_1{k+1}(3)*sin(xk_1{k+1}(4))*dt;
         dv_k_sol*dt;
         dtheta_k_sol*dt];
    f_xk_uk = f + xk_1{k+1};
     
    df_xk_uk = [1 0 cos(xk_1{k+1}(4))*dt  -xk_1{k+1}(3)*sin(xk_1{k+1}(4))*dt;
                0 1 sin(xk_1{k+1}(4))*dt  xk_1{k+1}(3)*cos(xk_1{k+1}(4))*dt;
                0 0 1 0; 
                0  0 0 1];
end

%Forward Pass

xk_opt{1} = x0; %Optimal initial state
for k = 1:N
    uk_opt = subs(u{k},[v_k theta_k dv_k dtheta_k ],[xk_opt{k}(3) xk_opt{k}(4) u{k}(1) u{k}(2)]); %Optimal Control Trajectory
    xk_opt{k+1} = xk_1{k};
    
end