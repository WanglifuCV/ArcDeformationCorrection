function S=constraints_Angles_curve(tau,X, Y,scale)
R=tau;
[m n]=size(X);
% U=X+(R.*ones(m,n)-(R.^2-Y.^2).^0.5);
% V=R.*asin(Y./R);
U=X+R.*(ones(m,n)-cos(Y./R));
V=R.*sin(Y./R);
e1 = [U(1)-U(3); V(1)-V(3)];
e2 = [U(2)-U(4); V(2)-V(4)];
norm_e1_2 = e1'*e1;
norm_e2_2 = e2'*e2;
e1e2 = e1'*e2;
N = 4 * sqrt( norm_e1_2 * norm_e2_2 - e1e2^2 );
dUdR=ones(m,n)-cos(Y./R)-Y./R.*sin(Y./R);
dVdR=-Y/R.*cos(Y/R);
S=constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N,dUdR,dVdR);