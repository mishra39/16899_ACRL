clear all 
clc
close all

Q = eye(2);
R = 1;
T = 0.1;

A = [1 T ;0 1];
B = [T^2 /2 ; T];

%% 3.2: Riccati Equation Solution
[K_dlqr,cov_dlqr, poles] = dlqr(A,B,Q,R);
[cov_dare,K_dare,L,info] = idare(A,B,Q,R,[],[]);

%% 3.4: Stochastic System Plot
sigma = 0.001;
sig_W = sqrt(0.0001);
sig_V = sqrt(0.00001);
k = 100;
x0 = [2;0];
x_k = [];
C = eye(2);
y_k_all = [];
for ind = 1:k
    w_k = normrnd(0,sqrt(sig_W));
    v_k = normrnd(0,sqrt(sig_V));
    u_k = -K_dlqr*x0;
    x_k_1 = A*x0 + B*(u_k + w_k);
    x_k = [x_k ; (x_k_1)'+v_k]; 
    x0 = x_k_1;
    y_k  = C*x_k_1 + v_k;
    y_k_all = [y_k_all ; y_k'];
end
figure (1)
plot(1:k,x_k)
title('Trajectory of the stochastic system')
xlabel('Time Steps')
ylabel('States')

%% 3.6: Steady State Kalman Gain and covariance
C = eye(2);
W = (0.0001);
V = 0.00001*eye(2);
Plant = ss(A,[B B], C, 0,-1);
[kalmf, L, P, M] = kalman(Plant,Q,R);
% [K_kf_dlqr,cov_kf_dlqr] = dlqr(A',C',W,V);
% [cov__kf_dare,K__kf_dare,L_kf,info_kf] = idare(A',C',W,V,[],[]);
[P_k, K_kf] = idare(A',C',Q,R);
%% 3.7
x0 = [2;0];
x_hat_k(:,1) = [2;0];
y_hat_k(:,1) = C*x_hat_k(:,1);
K_kf = P_k*C'*inv(C*P_k*C' + V);
for k = 2:100
    
    u_k = -K_dlqr*x_hat_k(:,k-1);
    x_hat_k(:,k) = x_hat_k(:,k-1) + (A*x_hat_k(:,k-1)+ B*u_k ) + K_kf*(y_k_all(k,:)' - y_hat_k(:,k-1));
    y_hat_k(:,k) = C*x_hat_k(:,k); 
end

figure (2)
plot(1:k ,x_hat_k)
title('Trajectory of the stochastic system')
xlabel('Time Steps')
ylabel('Predicted State')