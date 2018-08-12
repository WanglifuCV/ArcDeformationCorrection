function derivative = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdP, dVdP)
dNorm_e1_2dP = 2*e1(1)*(dUdP(1)-dUdP(3)) + 2*e1(2)*(dVdP(1)-dVdP(3));
dNorm_e2_2dP = 2*e2(1)*(dUdP(2)-dUdP(4)) + 2*e2(2)*(dVdP(2)-dVdP(4));
dE1e2_2dP = 2*e1e2* (e2(1)*(dUdP(1)-dUdP(3)) + e1(1)*(dUdP(2)-dUdP(4)) + e2(2)*(dVdP(1)-dVdP(3)) + e1(2)*(dVdP(2)-dVdP(4)) );

derivative = ( norm_e2_2*dNorm_e1_2dP +  norm_e1_2*dNorm_e2_2dP - dE1e2_2dP) / N;