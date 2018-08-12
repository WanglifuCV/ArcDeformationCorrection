 
    function [I1] = Transform_angle(I0, tau)

%f = 1250;
%sf=(2376.730859642686000+2347.583628880417600)/2; %1.48
% sf=(2433.522096879261900+2421.494747721276800 )/2;   %2.08
%sf=(2449.776340582253900 + 2436.347498311849900)/2;   %2.68
%sf=(1739.044847066450500+1730.166873006305400)/2;
%sf=(2350.804588244516700+2324.108345954108200)/2;;
%sf=(3422.15559+3413.41586)/2;
% sf=12.5/0.0052;
sf=8/0.0052;
tx = tau(1);
ty = tau(2);
tz = tau(3);
thetaZ = tau(4);
thetaY = tau(5);
thetaX = tau(6);
Sx = sin(thetaX); Cx = cos(thetaX);
Sy = sin(thetaY); Cy = cos(thetaY);
Sz = sin(thetaZ); Cz = cos(thetaZ);
R11 = Cy*Cz;
R21 = Sx*Sy*Cz - Cx*Sz;
R31 = Cx*Sy*Cz + Sx*Sz;
R12 = Cy*Sz;
R22 = Sx*Sy*Sz + Cx*Cz;
R32 = Cx*Sy*Sz - Sx*Cz;
t1 = tx * sf ;
t2 = ty * sf;
t3 = tz * sf;

[h, w] = size(I0);
cx = (w+1)/2;
cy = (h+1)/2;

I1 = zeros(h, w);
du=zeros(h,w);
dv=zeros(h,w);
for x = 1:w
    for y = 1:h
              
        xx = x - cx;
        yy = y - cx;
        
        N1 = f * R11 - R31 * xx;
        N2 = f^2 * tx - f * tz * xx;
        N3 = R32 * xx - f * R12;
        
        M1 = f * R22 - R32 * yy;
        M2 = f^2 * ty - f * tz * yy;
        M3 = R31 * yy - f * R21;
        
        
        xx = (N2*M1 + N3*M2) / (N3*M3 - N1*M1);
        yy = (N1*M2 + M3*N2) / (N3*M3 - N1*M1);
        
        xx = xx + cx;
        yy = yy + cy;
        
        x_1 = floor(xx);
        y_1 = floor(yy);
        x_2 = x_1 + 1;
        y_2 = y_1 + 1;
        x_s = min(max(x_1,1),w);
        y_s = min(max(y_1,1),h);
        f1 = I0(y_s,x_s);%(x_1,y_1)

        y_s = min(max(y_2,1),h);
        f2 = I0(y_s,x_s);%(x_1,y_2)

        x_s = min(max(x_2,1),w);
        y_s = min(max(y_1,1),h);
        f3 = I0(y_s,x_s);%(x_2,y_1)

        y_s = min(max(y_2,1),h);
        f4 = I0(y_s,x_s);%(x_2,y_2)


        f13 = (xx-x_1)*f3+(x_2-xx)*f1;
        f24 = (xx-x_1)*f4+(x_2-xx)*f2;

        I1(y, x) = (yy - y_1)*f24+(y_2-yy)*f13;
    end
end