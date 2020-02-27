function [t0, x0, u0] = P1_4_shift(T, t0, x0, u,f)
st = x0;
con = u(1,:)';
f_value = f(st,con);
st = st+ (T*f_value);
x0 = full(st);

t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)]; % Trim the first entry as we already used it in the previous iteration. Use last entry as a guess for next iteration
end