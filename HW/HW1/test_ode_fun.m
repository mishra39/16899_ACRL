function xdot=sys(t,X,K,A,B)
 u = K*x;
 xdot =A*X+B*U;
end