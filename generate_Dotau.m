function [Dotau, du, dv, X_loc_map, Y_loc_map] = generate_Dotau(img,tau, img_du, img_dv,XData,YData,image_center,scale)
% cx,cy：给定中心点，即变换的原点

[X0 Y0]=meshgrid(XData(1):XData(2), YData(1):YData(2));

[sx sy]=size(X0);


%sf = (2376.730859642686000+2347.583628880417600)/2; %1.48
%sf=(2433.522096879261900+2421.494747721276800 )/2;   %2.08
%sf=(2449.776340582253900 + 2436.347498311849900)/2;   %2.68
%sf=(1739.044847066450500+1730.166873006305400)/2;
%sf=(2350.804588244516700+2324.108345954108200)/2;;
%sf=(3422.15559+3413.41586)/2;
% sf=12.5/0.0052;
sf=8/0.0052;
sf=sf*scale;
tx = tau(1);
ty = tau(2);
tz = tau(3);
thetaZ = tau(4);
thetaY = tau(5);
thetaX = tau(6);
Sx = sin(thetaX); Cx = cos(thetaX);
Sy = sin(thetaY); Cy = cos(thetaY);
Sz = sin(thetaZ); Cz = cos(thetaZ);
R11 = Cy*Cz;  %
R21 = Sx*Sy*Cz - Cx*Sz;
R31 = Cx*Sy*Cz + Sx*Sz;
R12 = Cy*Sz;
R22 = Sx*Sy*Sz + Cx*Cz;
R32 = Cx*Sy*Sz - Sx*Cz;
t1 = tx * sf ;
t2 = ty * sf;
t3 = tz * sf;
% --------------------------------
[iHeight, iWidth]= size(img);

Dotau = zeros(sx,sy);
du = zeros(sx,sy);
dv = zeros(sx,sy);
X_loc_map = zeros(sx,sy);
Y_loc_map = zeros(sx,sy);

%K1 = R11*X0 + R12*Y0 + t1;   %% xx=img(i,j) for initial_img
%K2 = R21*X0  + R22*Y0 + t2;
%K3 = R31*X0  + R32*Y0 + t3;
%U2 = sf * K1/K3 ;       %%xx=img(i,j) for updated img  % cx,cy：给定中心点，即变换的原点
%V2 = sf * K2/K3 ;
U2=X0(1,:);
V2=Y0(:,1);

lu=length(U2);
lv=length(V2);

count=0;
%V2=V2+image_center(2)*ones(1,lv);
%U2=U2+image_center(1)*ones(lu,1);
for i=1:lu  %列数
    for j=1:lv  %行数
           
          
            K1 = R11*U2(i) + R12*V2(j) + t1;
            K2 = R21*U2(i) + R22*V2(j) + t2;
            K3 = R31*U2(i)+ R32*V2(j) + t3;
            xx = sf * K1/K3+image_center(1);
            yy = sf * K2/K3+image_center(2) ;
             if nargout > 1
                X_loc_map(j,i) =xx;%行
                Y_loc_map(j,i) =yy; %列
            end
            
            
            x_1 = floor(xx);
            y_1 = floor(yy);
            x_2 = x_1 + 1;
            y_2 = y_1 + 1;
            
            x_s = min(max(x_1,1),iWidth);
            y_s = min(max(y_1,1),iHeight);
            if x_s==1
                count=count+1;
            end
           % x_s=x_1;
            %y_s=y_1;
            f1 = img(y_s,x_s);%(x_1,y_1)
            if nargout > 1
                u1 = img_du(y_s,x_s);
                v1 = img_dv(y_s,x_s);
            end
            
            y_s = min(max(y_2,1),iHeight);
           % y_s=y_2;
            f2 = img(y_s,x_s) ; %(x_1,y_2)
            if nargout > 1
                u2 = img_du(y_s,x_s);
                v2 = img_dv(y_s,x_s);
            end
            
            x_s = min(max(x_2,1),iWidth);
            y_s = min(max(y_1,1),iHeight);
            %x_s=x_2;
            %y_s=y_1;
            f3 = img(y_s,x_s);  %(x_2,y_1)
            if nargout > 1
                u3 = img_du(y_s,x_s);
                v3 = img_dv(y_s,x_s);
            end
            
            y_s = min(max(y_2,1),iHeight);
            %y_s=y_2;
            f4 = img(y_s,x_s);   %(x_2,y_2)
            if nargout > 1
                u4 = img_du(y_s,x_s);
                v4 = img_dv(y_s,x_s);
            end
            
            a = (xx-x_1); b = (x_2-xx);
            c = (yy - y_1); d = (y_2-yy);
            
            f13 = a*f3+b*f1;
            f24 = a*f4+b*f2;
           
       
             Dotau(j,i) = c*f24+d*f13;
          
            
            if nargout > 1
                u13 = a*u3+b*u1;
                u24 = a*u4+b*u2;
                
                v13 = a*v3+b*v1;
                v24 = a*v4+b*v2;
                
                du(j,i) = c*u24+d*u13; 
                dv(j,i) = c*v24+d*v13;
            end
    end
end

    

