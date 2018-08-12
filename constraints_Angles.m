function S=constraints_Angles(tau,X, Y,scale)
% constraints() 约束变换前后变换区域的面积基本保持不变
% mode.
% -----------------------------input--------------------------------------
% X, Y:   4-by-1 real vector, respectively X,Y coordinates of 4 transform control
%                 points.
% tau:    变换参数
% ds:     像元尺寸
% scale:  当前处理的图像相对于原图像的缩放比例
% ----------------------------output--------------------------------------
% S:            linearized constraints on tau parameters.

% f =  tau(1);

%sf = (2376.730859642686000+2347.583628880417600)/2;  %1.48
%sf=(2433.522096879261900+2421.494747721276800 )/2;   %2.08
%sf=(2449.776340582253900 + 2436.347498311849900)/2;   %2.68
%sf=(1739.044847066450500+1730.166873006305400)/2;
%sf=(2350.804588244516700+2324.108345954108200)/2;;
%sf=(3422.15559+3413.41586)/2;
% sf=(1718.75056+1715.88691)/2;
% sf=12.5/0.0052;
sf=8/0.0052;
sf=sf*scale;
tx = tau(1);
ty = tau(2);
tz = tau(3);

Sx = sin(tau(6)); Cx = cos(tau(6));
Sy = sin(tau(5)); Cy = cos(tau(5));
Sz = sin(tau(4)); Cz = cos(tau(4));
R11 = Cy*Cz;
R21 = Sx*Sy*Cz - Cx*Sz;
R31 = Cx*Sy*Cz + Sx*Sz;
R12 = Cy*Sz;
R22 = Sx*Sy*Sz + Cx*Cz;
R32 = Cx*Sy*Sz - Sx*Cz;

K1 = R11*X + R12*Y + sf * tx;
K2 = R21*X + R22*Y + sf * ty;
K3 = R31*X + R32*Y + sf * tz;

% 变换后四个控制点坐标
U = sf * K1./K3;
V = sf * K2./K3;

e1 = [U(1)-U(3); V(1)-V(3)];
e2 = [U(2)-U(4); V(2)-V(4)];
norm_e1_2 = e1'*e1;
norm_e2_2 = e2'*e2;
e1e2 = e1'*e2;

K3_2 = K3.^2;
dK1dThetaZ = (-Cy*Sz)*X + (Cy*Cz)*Y;
dK1dThetaY = (-Sy*Cz)*X + (-Sy*Sz)*Y;
dK1dThetaX = 0;
dK2dThetaZ = (-Sx*Sy*Sz-Cx*Cz)*X + (Sx*Sy*Cz-Cx*Sz)*Y;
dK2dThetaY = (Sx*Cy*Cz)*X + (Sx*Cy*Sz)*Y;
dK2dThetaX = (Cx*Sy*Cz + Sx*Sz)*X + (Cx*Sy*Sz-Sx*Cz)*Y;
dK3dThetaZ = (Sx*Cz - Cx*Sy*Sz)*X + (Cx*Sy*Cz + Sx*Sz)*Y;
dK3dThetaY = (Cx*Cy*Cz)*X + (Cx*Cy*Sz)*Y;
dK3dThetaX = (Cx*Sz - Sx*Sy*Cz)*X + (-Sx*Sy*Sz-Cx*Cz)*Y;

N = 4 * sqrt( norm_e1_2 * norm_e2_2 - e1e2^2 );
S = zeros(1,6);

% % dS/df
% dUdf = s * (K1./K3 + sf * (tx * K3 - tz * K1)./K3_2);
% dVdf = s * (K2./K3 + sf * (ty * K3 - tz * K2)./K3_2);
% S(1, 1) = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdf, dVdf); 

% dS/dTx
dUdTx = (sf)^2./K3;
dVdTx = [0; 0; 0; 0];
S(1, 1) = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdTx, dVdTx);

% dS/dTy
dUdTy = dVdTx;
dVdTy = dUdTx;
S(1, 2) = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdTy, dVdTy);

% dS/dTz
dUdTz = -(sf)^2 * K1./K3_2;
dVdTz = -(sf)^2 * K2./K3_2;
S(1, 3) = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdTz, dVdTz);

% dS/dThetaZ
dUdThetaZ = sf * (dK1dThetaZ.*K3 - dK3dThetaZ.*K1)./K3_2;
dVdThetaZ = sf * (dK2dThetaZ.*K3 - dK3dThetaZ.*K2)./K3_2;
S(1, 4) = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdThetaZ, dVdThetaZ); 

% dS/dThetaY
dUdThetaY = sf * (dK1dThetaY.*K3 - dK3dThetaY.*K1)./K3_2;
dVdThetaY = sf * (dK2dThetaY.*K3 - dK3dThetaY.*K2)./K3_2;
S(1, 5) = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdThetaY, dVdThetaY); 

% dS/dThetaX
dUdThetaX = sf * (dK1dThetaX.*K3 - dK3dThetaX.*K1)./K3_2;
dVdThetaX = sf * (dK2dThetaX.*K3 - dK3dThetaX.*K2)./K3_2;
S(1, 6) = constrain_derivative(e1,e2, norm_e1_2, norm_e2_2, e1e2, N, dUdThetaX, dVdThetaX);