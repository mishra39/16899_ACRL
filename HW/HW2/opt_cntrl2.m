% function [u] = opt_cntrl2(x0,t,G_bar)
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

H = B_bar'*Q_bar*B_bar + R_bar;
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

x_exe = [];
x_pred = [];
x0 = [2;0];
for i = 1:N_total+50
    f_bar_x0 = f_bar*x0;
    ref_ind  = 0 ;% Slide window for reference
    f = (f_bar_x0'*Q_bar*B_bar) - (G_bar(2*(i)-1:2*(i-1)+2*(N_total+1),1)'*Q_bar*B_bar);
    A_quad = [L_bar_u; L_bar_xl*B_bar;L_bar_xu*B_bar];
    b_quad = [b_bar_u ; b_bar_xl - L_bar_xl*f_bar_x0; b_bar_xu - L_bar_xu*f_bar_x0];
    u = quadprog(H,f,A_quad,b_quad);
    u_exe = u(1); % Executed MPC
    x_k_1 = A*x0 + B*u_exe;
    x_exe = [x_exe ; x_k_1'];
    
    x_pred_i = f_bar_x0 + B_bar * u; % Predicted Trajectory
    x_pred = [x_pred ; x_pred_i'];
    x0 = x_k_1;
end



function [f_bar] = f_bar_mat(A,N_total)
        f_bar = [];
for i = 0:N_total
    f_int  = A^i;
    f_bar = cat(1,f_bar,f_int);
end
% f_bar(end-1,:) = 0
% f_bar(end,:) = 0
% disp(f_bar)
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
%   PHI(end-1,:) = 0
% PHI(end,:) = 0
end