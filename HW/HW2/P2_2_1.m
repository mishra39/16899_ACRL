%%
clear all
clc
close all
T = 0.1;
N = 10;
N_total = N/T;
A = [1 T; 0 1];
B = [(T^2)/2; T];
x0 = [2;0]; % Initial condition
G_ref = ref_st(N_total,N,x0);
x_pred = [];
x_act = [];
x_act = [x_act; x0'];
x_pred_all = [];
for k = 0:N_total
    [u,fval,exitflag,output,lambda] = opt_cntrl2(x0,k,G_ref);
    x_k_1 = A*x0 + B*(u(1)+ mvnrnd(0,0.001)); % Executed Trajectory
    x_act =  [x_act ; x_k_1'];
    
    % Predicted Trajectory
    for m = 1:size(u,1)
        x_pred_val = A*x0 + B*(u(m) + mvnrnd(0,0.001));
        x_pred = [x_pred;x_pred_val'];
    end
    x_pred_all = [x_pred_all ; x_pred];
    x_pred = [];
    x0= x_k_1;
end
% Predicted
pred_x1 = x_pred_all(:,1);
pred_x2 = x_pred_all(:,2);
% Executed
act_x1 = x_act(:,1);
act_x2 = x_act(:,2);

% Reference
ref_x1 = [];
ref_x2 = [];
ref_x1 = [ref_x1;G_ref(1,1)]

for i = 1:(size((G_ref),1)/2)-1
    ref_x11 = G_ref((2*i)+1,1);
    ref_x22 = G_ref(2*i,1);
    ref_x1 = [ref_x1;ref_x11];
    ref_x2 = [ref_x2;ref_x22];
end
ref_x2 = [ref_x2;G_ref(end)];
%%
figure(1)
% plot(pred_x1,pred_x2,'--')
% hold on
% plot(act_x1,act_x2,'o')
% hold on
plot(ref_x1,ref_x2,'x')
xlabel('x1')
ylabel('x2')
ylim([-2 2])
title('Noisy Predicted and Executed Trajectories')
legend('Predicted','Executed','Reference')