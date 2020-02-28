function [x_pred_all,x_act] = MPC_calc(w_k)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
T = 0.1;
N = 10;
N_total = N/T;
A = [1 T; 0 1];
B = [(T^2)/2; T];
x0 = [2;0]; % Initial condition
G_ref = ref_st(N_total,N,x0)
x_pred = [];
x_act = [];
x_act = [x_act; x0'];
x_pred_all = [];

for k = 0:N_total+100
    u = opt_cntrl2(x0,k,G_ref);
    if w_k == 0
         x_k_1 = A*x0 + B*(u(1));
    else
        x_k_1 = A*x0 + B*(u(1)+ mvnrnd(0,0.001)); % Executed Trajectory
    end
    x_act =  [x_act ; x_k_1'];
    
    % Predicted Trajectory
    for m = 1:size(u,1)
        if w_k ==0
            x_pred_val = A*x0 + B*(u(m));
        else
            x_pred_val = A*x0 + B*(u(m) + mvnrnd(0,0.001));
        end
        
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
ref_x1 = [ref_x1;G_ref(1,1)];

for i = 1:(size((G_ref),1)/2)-1
    ref_x11 = G_ref((2*i)+1,1);
    ref_x22 = G_ref(2*i,1);
    ref_x1 = [ref_x1;ref_x11];
    ref_x2 = [ref_x2;ref_x22];
end
disp(ref_x1)
ref_x2 = [ref_x2;G_ref(end)];
ref_x1_plot_all = [ref_x1;zeros(size(act_x1,1)-size(ref_x1,1),1)];
%%
figure(1)
% plot(pred_x1, pred_x2,'--')
hold on
plot(1:size(act_x1),act_x1,'o')
plot(1:size(ref_x1),ref_x1,'--')
% plot(1:sizeact_x1,'o')

% hold on
% plot(1:size(ref_x1),G_ref(,'x')
xlabel('t')
ylabel('actual_x1')
title('Noisy Predicted and Executed Trajectories')

legend('Executed','Reference')
end