function [Dotau] = generate_new_image(img,tau,image_center,scale)

%tau= load('tau.mat');
%scale=load('scale.mat');
%[img_name, img_path]=uigetfile('*.*');
%img=imread(fullfile(img_path, img_name));
%[H W]=size(img); %H 是行数 
[iHeight, iWidth]= size(img);


[X0 Y0]=meshgrid(1:iWidth, 1:iHeight); %X是行坐标
[sx sy]=size(X0);
U2=X0(1,:);
V2=Y0(:,1);
U2=U2-iWidth/2*ones(1,iWidth);
V2=V2-iHeight/2*ones(iHeight,1);
lu=length(U2);
lv=length(V2);
%scale=scale.scale
%sf = 2.4375e+003;
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
Dotau = zeros(sx,sy);

for i=1:lu  %列数
    for j=1:lv  %行数  
            K1 = R11*U2(i) + R12*V2(j) + t1;
            K2 = R21*U2(i) + R22*V2(j) + t2;
            K3 = R31*U2(i)+ R32*V2(j) + t3;
            xx = sf * K1/K3+image_center(1);
            yy = sf * K2/K3+image_center(2) ;   
            
            x_1 = floor(xx);
            y_1 = floor(yy);
            x_2 = x_1 + 1;
            y_2 = y_1 + 1;
            
           if(x_1<1||y_1<1||x_2>iWidth||y_2>iHeight)
               Dotau(i,j)=0;
           else           
            x_s = min(max(x_1,1),iWidth);
            y_s = min(max(y_1,1),iHeight);
            f1 = img(y_s,x_s);%(x_1,y_1)
   
            
            y_s = min(max(y_2,1),iHeight);
            f2 = img(y_s,x_s); %(x_1,y_2)

            
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
end
%figure;imshow(uint8(Dotau));