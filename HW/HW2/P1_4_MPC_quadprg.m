clc
clear all
close all
T = 0.1;
N = 10;
A = [1 T; 0 1];
B = [(T^2)/2; T];
C = [1 0];
D = 0;
Q = eye(2);
R = 1;
S = 10;
f_bar = [];
B_bar = [];
N_total = N/T;

%%
Q_bar = eye((N_total+1)*2);
Q_bar(end) = S;
Q_bar(end-1,end-1) = S;

R_bar = eye(N_total);
G_bar = [];
x0 = [2;0];
for k = 1:(N_total+1)*3
    if k<=N_total
        ref = ((N_total-(k-1))/N_total)*x0;
        G_bar = cat(1,G_bar,ref);

    else
        ref = [0 ;0];
        G_bar = cat(1,G_bar,ref);
    end
end

f_bar = f_bar_mat(A,N_total);
B_bar = B_bar_matrix(A,B,N_total,N_total);

H = 0.5* B_bar'*Q_bar*B_bar + R_bar;
L_bar_ul = -1*eye(N_total);
L_bar_uu = eye(N_total);
L_bar_u = cat(1,L_bar_ul,L_bar_uu);
b_bar_ul = ones(N_total,1);
b_bar_uu = ones(N_total,1);
b_bar_u = cat(1,b_bar_ul,b_bar_uu);

L_bar_xl = -1*eye((N_total+1)*2);
L_bar_xu = eye((N_total+1)*2);
L_bar_x = cat(1,L_bar_xl,L_bar_xu);
b_bar_xl = 5*ones((N_total+1)*2,1);
b_bar_xu = 5*ones((N_total+1)*2,1);
b_bar_x = cat(1,b_bar_xl,b_bar_xu);

%% 1.2: Feasibility Map 
x1_limit = -5:0.1:5;
x2_limit = -5:0.1:5;
[X,Y] = meshgrid(x1_limit, x2_limit);
feas_mat = []; % feasibilty matrix

for st1 = 1:size(X,1)
    st22 = 0 ;
    for st2 = size(Y,1):-1:1
        st22  = st22 + 1;
        x1 = X(st2,st1);
        x2 = Y(st2,st1);
        x2_0 = [x1;x2]
        f_bar2_x0 = f_bar*x2_0;
        f2 = (f_bar2_x0'*Q_bar*B_bar) - (G_bar(2*(1)-1:2*(1-1)+2*(N_total+1),1)'*Q_bar*B_bar);
        A_quad2 = [L_bar_u; L_bar_xl*B_bar;L_bar_xu*B_bar];
        b_quad2 = [b_bar_u ; b_bar_xl - L_bar_xl*f_bar2_x0; b_bar_xu - L_bar_xu*f_bar2_x0];
        [u2,fval,exitflag,output,lambda] =  quadprog(H,f2,A_quad2,b_quad2);
        if exitflag > 0
            feas_mat(st22,st1) = 1;
        else
            feas_mat(st22,st1) = 0;
        end
    end
end
%% Saving Parameters
figure(1)
plt_map = image(feas_mat*255);
colormap
xlim([-5 5])
filename = 'feasibility_final.mat';
save(filename,'feas_mat');

%% 1.4
x_exe = [];
x_pred = [];
u_opt = [];
x0 = [2;0];
for i = 1:N_total+50
    f_bar_x0 = f_bar*x0;
    f = (f_bar_x0'*Q_bar*B_bar) - (G_bar(2*(i)-1:2*(i-1)+2*(N_total+1),1)'*Q_bar*B_bar);
    A_quad = [L_bar_u; L_bar_xl*B_bar;L_bar_xu*B_bar];
    b_quad = [b_bar_u ; b_bar_xl - L_bar_xl*f_bar_x0; b_bar_xu - L_bar_xu*f_bar_x0];
    u = quadprog(H,f,A_quad,b_quad);
    u_exe = u(1); % Executed MPC
    u_opt = [u_opt ; u_exe];
    x_k_1 = A*x0 + B*u_exe;
    x_exe = [x_exe ; x_k_1'];
    
    x_pred_i = f_bar_x0 + B_bar * u; % Predicted Trajectory
    x_pred = [x_pred ; x_pred_i'];
    figure(2)
%     plot(x_pred_i(1:2:end), x_pred_i(2:2:end))
%     hold on
    x0 = x_k_1;
end
%Plotting
plot(x_exe(:,1), x_exe(:,2))
xlim([-0.5 2.2])
xlabel('Number of steps')
ylabel('Value of parameter')
title('Executed and Predicted Trajectories')
legend('Predicted','Executed')

%%
figure(3)
plot(1:size(x_exe,1), x_pred,'b-')
hold on
plot(1:size(x_exe,2),x_exe,'go')
xlabel('Number of steps')
ylabel('Value of parameter')
title('Executed and Predicted Trajectories')
legend('Predicted','Executed')
%% 2.1: MPC with noise

x_exe2 = [];
x_pred2 = [];
x0 = [2;0];
sigma = sqrt(0.00001);
mu = 0;
for i = 1:N_total+50
    f_bar2_x0 = f_bar*x0;
    w_k = normrnd(mu,sigma);
    f2 = (f_bar2_x0'*Q_bar*B_bar) - (G_bar(2*(i)-1:2*(i-1)+2*(N_total+1),1)'*Q_bar*B_bar);
    A_quad2 = [L_bar_u; L_bar_xl*B_bar;L_bar_xu*B_bar];
    b_quad2 = [b_bar_u ; b_bar_xl - L_bar_xl*f_bar2_x0; b_bar_xu - L_bar_xu*f_bar2_x0];
    u = quadprog(H,f2,A_quad2,b_quad2);
    u_exe = u(1); % Executed MPC
    x_k_1 = A*x0 + B*(u_exe+w_k);
    x_exe2 = [x_exe2 ; x_k_1'];
    
    x_pred_i = f_bar_x0 + B_bar * (u+w_k); % Predicted Trajectory
    x_pred2 = [x_pred2 ; x_pred_i'];
    x0 = x_k_1;
    figure(4)
    plot(x_pred_i(1:2:end), x_pred_i(2:2:end))
    hold on
end

%% Plotting
plot(x_exe2(:,1), x_exe2(:,2))
xlim([-2 2])
xlabel('Number of steps')
ylabel('Value of parameter')
title('Executed and Predicted Trajectories with Noise')
legend('Predicted','Executed')

figure(5)
% plot(1:size(x_exe2,1), x_pred2,'b-')
% hold on
plot(1:size(x_exe2,1),x_exe2,'r-')
xlabel('Number of steps')
ylabel('Value of parameter')
title('Executed Trajectory with Noise')
legend('Predicted','Executed')
%% 2.2: Error from trajectory
E_1 = x_exe - x_exe2;

%% 2.3: Design L_k gain
uk_mpc = u_opt;
x1_opt = x_exe;
x0 = [2;0];
L_k0  = zeros(1,2);
L_k_all = [];
opt_xk = [];
f_val_mat = [];
Q = 0.9;
for i = 1:size(uk_mpc,1)
    error_2 = (@(L) norm((x1_opt(i,:) - (A*x0 + Q*B*(uk_mpc(i) + L(1)*E_1(i,1)+ L(2)*E_1(i,2))))));
    [L,f_val] = fmincon(error_2,L_k0,[],[]);
    x0 = (A*x0 + B*(uk_mpc(i) + L(:,1)*E_1(i,1)+ L(:,2)*E_1(i,2)));
    L_k_all = [L_k_all; L];
    f_val_mat  = [f_val_mat ; f_val];
end
figure(6)
plot(f_val_mat,'*')
xlabel('Number of steps')
ylabel('Error Magnitude')
title('Error with the new Control Law')


%% 2.4: Feedforward Learning Gain
x_exe2 = [];
x_pred2 = [];
x0 = [2;0];
sigma = 0.007;
mu = 0;
Q = 0.9;
x_exec_L_Q_mat = [];
for i = 1:N_total+50
    x_exec_L_Q =  A*x0 + B*(Q*(uk_mpc(i) + L_k_all(i,1)*E_1(i,1)+ L_k_all(i,2)*E_1(i,2)));
    x0 = x_exec_L_Q;
    x_exec_L_Q_mat = [x_exec_L_Q_mat; x_exec_L_Q'];
end
figure(7)
plot(x_exec_L_Q_mat,'b-')
xlabel('Number of steps')
ylabel('Value of parameter')
title('Executed Trajectory with Feedforward Gain and Q-filter')

%%
function out = E_2(L,x1_opt, A, B, uk_mpc, E_1,x0,i)
    disp(L)
    error_2 = (x1_opt(i,:) - (A*x0 + B*(uk_mpc(i) + L(1)*E_1(i,1)+ L(2)*E_1(i,2))));
    out = norm(error_2);
end
function [f_bar] = f_bar_mat(A,N_total)
        f_bar = [];
for i = 0:N_total
    f_int  = A^i;
    f_bar = cat(1,f_bar,f_int);
end
end

function [PHI] = B_bar_matrix(A, B, Np, Nc)
  F = [];
  PHI = [];
  for j = 1:Nc
    for i = (1-j):(Np-j)
      if i < 0
        F = [F; 0*A^i*B];
      else
        F = [F; A^i*B];
      end
    end
    % Add to PHI
    PHI = [PHI F];
    % Clear F
    F = [];
  end
  PHI = cat(1,zeros(2,Np),PHI);
end