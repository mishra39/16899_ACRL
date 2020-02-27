function [Q] = Q_mat(Q_val,N_total)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
F = [];
Q = [];
for i = 1:N_total
%     disp(i)
    for j = 1:N_total
        if (i==j)
            F = [F;Q_val];
        else
            F = [F;0*Q];
        end
        Q = [Q;F];
        F = [];
    end
end
end
