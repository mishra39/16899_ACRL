clear all 
clc
close all

Q = eye(2);
R = 1;
T = 0.1;

A = [1 T ;0 1];
B = [T^2 /2 ; T];

[P,K,L,info] = idare(A,B,Q,R,[],[]);
[K,S,e] = dlqr(A,B,Q,R);

sigma = 0.001;
sig_W = 0.0001;
sig_V = 0.00001;
w_k = mvnrnd(0,sig_W);
v_k = mvnrnd(0,sig_V);
k = 100;
x0 = [2;0];
x_k = [];

for ind = 1:k
    w_k = mvnrnd(0,sig_W);
    v_k = mvnrnd(0,sig_V);
    x_k = [x_k ; x0'];
    u_k = -K*x0;
%     y_k(k) = x_k(k) + v_k;
    x_k_1 = A*x0 + B*(u_k + w_k);
    x0 = x_k_1;
end
figure (1)
plot(1:k,x_k)
title('Trajectory of the stochastic system')
xlabel('Time Steps')
ylabel('Predicted State')

%% 3.6
C = eye(2);
V = 0.00001*eye(2);
E = B;
W = 0.0001;
B_w = C;
V = 0.00001*eye(2);
[M_s,KF_gain,L,info] = idare(A',C',B_w*W*B_w',[],V,[],[]);
Plant = ss(A,[B B],C,0,-1,'inputname',{'u' 'w'},'outputname','y');
[kalmf,L,P,M] = kalman(Plant,W,V);
kalmf = kalmf(1,:);
%% 3.7
A = [1 T ;0 1];
B = [T^2 /2 ; T];
[X,K,L,info] = idare(A,B,Q,R,[],[]);
W = 0.0001;
V = 0.00001*eye(2);
Bw_k_1 = C;
sys = ss(A,B,C);
x0 = [2;0];
x_k_1_k_1 = x0;
u_k_1 = 0;
M_k_minus_1 = 0;

for k = 1:20
    x_k_k_1 = A * x_k_1_k_1 + B * u_k_1;
    M_k = A * M_k_minus_1 * A' + Bw_k_1 * W * Bw_k_1'; % B_w???????
    F_k = M_k * C' * inv( C * M_k * C' + V);
    M_k_plus_1 =  M_k - F_k * C * M_k;
    x_k_k = x_k_k_1 + F_k *(y_k - (C_k * x_k_k_1)); % y_k???????
    u_k = -K * x_k_k;
    x_k_k_1 = x_k_k;
    M_k = M_k_plus_1;

end