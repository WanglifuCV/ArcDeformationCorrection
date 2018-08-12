
tau= load('tau.mat');
scale=load('scale.mat');
[img_name, img_path]=uigetfile('*.*');
img=imread(fullfile(img_path, img_name));
%[H W]=size(img); %H 是行数 
[iHeight, iWidth]= size(img);


[X0 Y0]=meshgrid(1:iWidth, 1:iHeight); %X是行坐标
[sx sy]=size(X0);
U2=X0(1,:);
V2=Y0(:,1);
U2=U2-250*ones(1,500);
V2=V2-250*ones(500,1);
lu=length(U2);
lv=length(V2);
scale=scale.scale
sf =12.5/0.0052;
sf=sf*scale;
tx = tau.tau(1);
ty = tau.tau(2);
tz = tau.tau(3);
thetaZ = tau.tau(4);
thetaY = tau.tau(5);
thetaX = tau.tau(6);
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


Dotau = zeros(sx,sy);
for i=1:lu  %列数
    for j=1:lv  %行数  
            K1 = R11*U2(i) + R12*V2(j) + t1;
            K2 = R21*U2(i) + R22*V2(j) + t2;
            K3 = R31*U2(i)+ R32*V2(j) + t3;
            xx = sf * K1/K3+250;             %+image_center(1);
            yy = sf * K2/K3+250;             %+image_center(2) ;   
            
            x_1 = floor(xx);
            y_1 = floor(yy);
            x_2 = x_1 + 1;
            y_2 = y_1 + 1;
            
            x_s = min(max(x_1,1),iWidth);
            y_s = min(max(y_1,1),iHeight);
            f1 = img(y_s,x_s);%(x_1,y_1)
   
            
            y_s = min(max(y_2,1),iHeight);
            f2 = img(y_s,x_s) ; %(x_1,y_2)

            
            x_s = min(max(x_2,1),iWidth);
            y_s = min(max(y_1,1),iHeight);
            f3 = img(y_s,x_s);  %(x_2,y_1)
            
            
            y_s = min(max(y_2,1),iHeight);
            f4 = img(y_s,x_s);   %(x_2,y_2)

            
            a = (xx-x_1); b = (x_2-xx);
            c = (yy - y_1); d = (y_2-yy);
            
            f13 = a*f3+b*f1;
            f24 = a*f4+b*f2;
           
       
             Dotau(j,i) = c*f24+d*f13;
          
    end
end
figure;imshow(uint8(Dotau));