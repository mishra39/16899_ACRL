function dxdt = func(t,x,u)

dxdt(1) = x(3)*cos(x(4));
dxdt(2) = x(3)*sin(x(4));
dxdt(3) = u(1);
dxdt(4) = u(2);

dxdt = dxdt(:);
end