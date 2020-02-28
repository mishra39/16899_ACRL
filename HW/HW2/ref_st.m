function [G_bar] = ref_st(N_total,N,x0)
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
disp(G_bar)
end

