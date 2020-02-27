% point stabilization + Single shooting
clear all
close all
clc

% CasADi v3.4.5
addpath('C:\Users\sunny\OneDrive\Documents\casadi-windows-matlabR2016a-v3.5.1')
import casadi.*

T = 0.1; %[s]
N = 10; % prediction horizon
rob_diam = 0.3;

v_max = 1; v_min = -v_max;
x = SX.sym('x'); y = SX.sym('y'); %theta = SX.sym('theta');
states = [x;y]; n_states = length(states);

v = SX.sym('v'); %omega = SX.sym('omega');
controls = [v]; n_controls = length(controls);
rhs = [x + (T*y) + ((T^2/2)*v); y + (T*v)]; % system r.h.s


f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);
% parameters (which include the initial and the reference state of the robot)

X = SX.sym('X',n_states,(N+1));
% A Matrix that represents the states over the optimization problem.
% compute solution symbolically
X(:,1) = P(1:2); % initial state

% Solving for system states over the horizon using the discretized system
% model
for k = 1:N
    st = X(:,k);  con = U(:,k);
    f_value  = f(st,con);
    st_next  = f_value; % Symbolic Value of the next state
    X(:,k+1) = st_next; % Forward Propagation of the state
end

% this function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,P},{X});

obj = 0; % Objective function
g = [];  % constraints vector

Q = eye(2);%,2); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.1; % weighing matrices (states)
R = 1; %zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighing matrices (controls)
S = 10*eye(2);
% compute objective
for k=1:N-1
    st = X(:,k);  con = U(:,k);
    obj = obj+(1/2)*(st-P(3:4))'*Q*(st-P(3:4)) + (1/2)*con'*R*con; % calculate obj in symbolic form of the optimal state
end
st = X(:,N);
obj = obj+(1/2)*(st-P(3:4))'*S*(st-P(3:4));
% compute constraints
for k = 1:N+1   % box constraints to stay within the map margins
    g = [g ; X(1,k)];   %state x
    g = [g ; X(2,k)];   %state y
end

% make the decision variables one column vector
OPT_variables = reshape(U,N,1); % Two decision variables at each step of the horizon: inputs v and omega
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P); % Defining the Non-linear programming problem. 
% ******************x is now the vector of decision variables, i.e. control inputs***********************

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
% inequality constraints (state constraints)
args.lbg = -5;  % lower bound of the states x and y
args.ubg = 5;   % upper bound of the states x and y 

% input constraints
args.lbx(1:N,1) = v_min; %args.lbx(2:2:2*N,1)   = omega_min;
args.ubx(1:N,1) = v_max; %args.ubx(2:2:2*N,1)   = omega_max;
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [2 ; 0];    % initial condition.
%  xs = [1.5 ; 1.5]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,1);  % two control inputs 

sim_tim = 20; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 1e-2 && mpciter < sim_tim / T)
    if(mpciter <= N)
        xs = [((N-mpciter)/N)*x0(1) ; ((N-mpciter)/N)*x0(2)]; % Reference posture
        
    else
        xs = [0 ; 0];
        
    end
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = reshape(u0',N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)',1,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:2,mpciter+1)= full(ff_value)';
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0] = P1_4_shift(T, t0, x0, u,f); % get the initialization of the next optimization step
    xx(:,mpciter+2) = x0;  
    mpciter = mpciter + 1;
end;
main_loop_time = toc(main_loop);
ss_error = norm((x0-xs),2)
average_mpc_time = main_loop_time/(mpciter+1)
%%
% Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam,'Single_Shooting_N_150') % a drawing function