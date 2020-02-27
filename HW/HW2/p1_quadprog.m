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

[u,fval,exitflag,output,lambda] = opt_cntrl(x0);
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
x0 = [2;0];
x_pred = [];
x_act = [];
x_act = [x_act; x0'];
x_pred_all = [];
for k = 0:N_total
    [u,fval,exitflag,output,lambda] = opt_cntrl2(x0,k);
    x_k_1 = A*x0 + B*u(1); % Executed Trajectory
    x_act =  [x_act ; x_k_1'];
    
    % Predicted Trajectory
    for m = 1:size(u,1)
        x_pred_val = A*x0 + B*u(m);
        x_pred = [x_pred;x_pred_val'];
    end
    x_pred_all = [x_pred_all ; x_pred];
    x_pred = [];
    x0= x_k_1;
end