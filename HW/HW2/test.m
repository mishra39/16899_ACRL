F = [];
Q = [];
Q_val = eye(2);
N_total = 4;
for i = 1:N_total
%     disp(i)
    for j = (1-i):(N_total-i)
        if (i==j)
            F = [F;Q_val];
        else
            F = [F;0*Q];
        end
        Q = [Q;F];
        F = [];
    end
end