function [u,fval,exitflag,output,lambda] = opt_cntrl2(x0,t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
function [f_bar] = f_bar_mat(A,N_total,x0)
        f_bar = [];
for i = 0:N_total
    f_int  = A^i;
    f_bar = cat(1,f_bar,f_int);
end
f_bar = f_bar *x0;
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
Q_bar = eye((N_total+1)*2);
Q_bar(end) = S;
Q_bar(end-1,end-1) = S;

R_bar = eye(N_total);

G_bar = [];%zeros(size(f_bar));
for k = 1:(N_total+1)*3
    if k<=N
        ref = ((N-k)/N)*x0;%[0;0]
        G_bar = cat(1,G_bar,ref);
    else
        ref = [0 ;0];
        G_bar = cat(1,G_bar,ref);
    end
end
G_bar_use = G_bar(t+1:t+1++1+(N_total*2));

disp(size(G_bar_use))
f_bar = f_bar_mat(A,N_total,x0);
B_bar = B_bar_matrix(A,B,N_total,N_total);

H = B_bar'*Q_bar*B_bar + R_bar;
f = (f_bar'*Q_bar*B_bar) - (G_bar_use'*Q_bar*B_bar);%+ (1/2)*B_bar'*Q_bar*B_bar ;
disp(size(H))
disp(size(f))
disp(size(f_bar))
disp(size(B_bar))
disp(size(Q_bar))
disp(size(R_bar))
disp(size(G_bar_use))
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

% Constraints
A_quad = [L_bar_u; L_bar_xl*B_bar;L_bar_xu*B_bar];
b_quad = [b_bar_u ; b_bar_xl - L_bar_xl*f_bar; b_bar_xu - L_bar_xu*f_bar];

% A_quad = cat(1, L_bar_ul, L_bar_uu,L_bar_xl, L_bar_xu);
% b_quad = cat(1, b_bar_ul, b_bar_uu,b_bar_xl, b_bar_xu);
[u,fval,exitflag,output,lambda] = quadprog(H,f,A_quad,b_quad);
end