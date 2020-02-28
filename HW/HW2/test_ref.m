T = 0.1;
N = 10;
N_total = N/T;
A = [1 T; 0 1];
B = [(T^2)/2; T];
x0 = [2;0]; % Initial condition
for k = 1:(N_total+1)*3
    if k<=N_total
        ref = ((N_total-(k-1))/N_total)*x0;%[0;0]
        G_bar = cat(1,G_bar,ref);

    else
        ref = [0 ;0];
        G_bar = cat(1,G_bar,ref);
    end
end

x_pred = [];
x_act = [];
x_act = [x_act; x0'];
x_pred_all = [];