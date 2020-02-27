function [u,fval,exitflag,output,lambda] = opt_cntrl(x0)
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

for i = 0:N
    f_int  = A^i;
    f_bar = cat(1,f_bar,f_int);
end

f_bar = f_bar *x0
G_bar = [];%zeros(size(f_bar));
for k = 0:N
    ref = ((N-k)/N)*x0;%[0;0]
    G_bar = cat(1,G_bar,ref);
        
end


B_bar = [zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1)  zeros(2,1);
    B zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1)  zeros(2,1);
        A*B B zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1);
        A^2*B A*B B zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1);
        A^3*B A^2*B A*B B  zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1);
        A^4*B A^3*B A^2*B A*B B zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1);
        A^5*B A^4*B A^3*B A^2*B A*B B zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1);
        A^6*B A^5*B A^4*B A^3*B A^2*B A*B B zeros(2,1) zeros(2,1) zeros(2,1);
        A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B A*B B zeros(2,1) zeros(2,1) ;
        A^8*B A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B A*B B zeros(2,1);
        A^9*B A^8*B A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B A*B B;];
    
Q_bar = [Q zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2)  zeros(2,2) zeros(2,2);
    zeros(2,2) Q zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
    zeros(2,2) zeros(2,2) Q zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
    zeros(2,2) zeros(2,2) zeros(2,2) Q zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);
    zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) Q zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) ;
    zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) Q zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) ;
    zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) Q zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) ;
    zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) Q zeros(2,2) zeros(2,2) zeros(2,2) ;
    zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) Q zeros(2,2) zeros(2,2) ;
    zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) Q  zeros(2,2);
    zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) Q ];

R_bar = eye(10);

H = B_bar'*Q_bar*B_bar + R_bar;
f = (f_bar'*Q_bar*B_bar) - (G_bar'*Q_bar*B_bar);%+ (1/2)*B_bar'*Q_bar*B_bar ;
% disp(size(f))
L_bar_ul = -1*eye(10);
L_bar_uu = eye(10);
L_bar_u = cat(1,L_bar_ul,L_bar_uu);
b_bar_ul = ones(10,1);
b_bar_uu = ones(10,1);
b_bar_u = cat(1,b_bar_ul,b_bar_uu);


L_bar_xl = -1*eye(22);
L_bar_xu = eye(22);
L_bar_x = cat(1,L_bar_xl,L_bar_xu);
b_bar_xl = 5*ones(22,1);
b_bar_xu = 5*ones(22,1);
b_bar_x = cat(1,b_bar_xl,b_bar_xu);

% Constraints
A_quad = [L_bar_u; L_bar_xl*B_bar;L_bar_xu*B_bar];
b_quad = [b_bar_u ; b_bar_xl - L_bar_xl*f_bar; b_bar_xu - L_bar_xu*f_bar];

disp(size(H))
disp(size(f))
disp(size(f_bar))
disp(size(B_bar))
disp(size(Q_bar))
disp(size(R_bar))
disp(size(G_bar))
% A_quad = cat(1, L_bar_ul, L_bar_uu,L_bar_xl, L_bar_xu);
% b_quad = cat(1, b_bar_ul, b_bar_uu,b_bar_xl, b_bar_xu);
[u,fval,exitflag,output,lambda] = quadprog(H,f,A_quad,b_quad);
end