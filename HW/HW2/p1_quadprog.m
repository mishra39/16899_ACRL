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
% feas_mat = load(filename);

% [u,fval,exitflag,output,lambda] = opt_cntrl(x0);
x1_limit = -5:0.1:5;
x2_limit = -5:0.1:5;
[X,Y] = meshgrid(x1_limit, x2_limit);
feas_mat = []; % feasibilty matrix
%% 1.2: Feasibility Map 
for st1 = 1:size(X,1)
    for st2 = 1:size(Y,1)
        x1 = X(st2,st1);
        x2 = Y(st2,st1);
        x0 = [x1;x2]
        [u,fval,exitflag,output,lambda] = opt_cntrl2(x0,0);
        if exitflag == 1
            feas_mat(st2,st1) = 1;
        else
            feas_mat(st2,st1) = 0;
        end
    end
end
%% Saving Parameters
filename = 'feasibility.mat';
save(filename,'feas_mat');

%% 1.4
% Define reference
G_bar = [];
for k = 1:(N_total+1)*3
    if k<=N_total
        ref = ((N_total-(k-1))/N_total)*x0;%[0;0]
        G_bar = cat(1,G_bar,ref);

    else
        ref = [0 ;0];
        G_bar = cat(1,G_bar,ref);
    end
end

  F = [];
  PHI = [];
  for j = 1:N_total
    for i = (1-j):(N_total-j)
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
  PHI = cat(1,zeros(2,N_total),PHI);
  B_bar = [];
  B_bar = PHI;
%%




w=0;
[x_pred1,x_act1] = MPC_calc(w);
%%
w = 10;
[x_pred2,x_act2] = MPC_calc(w);

% Error
E_1 = x_act1 - x_act2;

% x0 = [2;0];
% G_ref = ref_st(N_total,N,x0);
% x_pred = [];
% x_act = [];
% x_act = [x_act; x0'];
% x_pred_all = [];
% for k = 0:N_total
%     [u,fval,exitflag,output,lambda] = opt_cntrl2(x0,k,G_ref);
%     x_k_1 = A*x0 + B*u(1); % Executed Trajectory
%     x_act =  [x_act ; x_k_1'];
%     
%     % Predicted Trajectory
%     for m = 1:size(u,1)
%         x_pred_val = A*x0 + B*u(m);
%         x_pred = [x_pred;x_pred_val'];
%     end
%     x_pred_all = [x_pred_all ; x_pred];
%     x_pred = [];
%     x0= x_k_1;
% end
% % Predicted
% pred_x1 = x_pred_all(:,1);
% pred_x2 = x_pred_all(:,2);
% % Executed
% act_x1 = x_act(:,1);
% act_x2 = x_act(:,2);
% 
% % Reference
% ref_x1 = [];
% ref_x2 = [];
% ref_x1 = [ref_x1;G_ref(1,1)];
% 
% for i = 1:(size((G_ref),1)/2)-1
%     ref_x11 = G_ref((2*i)+1,1);
%     ref_x22 = G_ref(2*i,1);
%     ref_x1 = [ref_x1;ref_x11];
%     ref_x2 = [ref_x2;ref_x22];
% end
% ref_x2 = [ref_x2;G_ref(end)];
% %%
% figure(1)
% plot(pred_x1,pred_x2,'--')
% hold on
% plot(act_x1,act_x2,'o')
% hold on
% plot(ref_x1,ref_x2,'x')
% xlabel('x1')
% ylabel('x2')
% ylim([0 2])
% title('Predicted and Executed Trajectories')
% legend('Predicted','Executed','Reference')