clear;
clc;
% y=-1:0.01:1;
% x=0;
% plot(x,y,'r-');
% R=3;
% x1=x-(R-(R^2-y.^2).^0.5);
% y1=R*asin(y./R);
% plot(x1,y1,'r-');

% R=3;
% y1=-1:0.01:1;
% y=R*sin(y1./R);
% x1=0;
% x=x1+(R-(R^2-y1.^2).^0.5);
% plot(x1,y1,'b-');hold on;
% plot(x,y,'r-');hold off;
% axis([-1 1 -1 1]);
% x2=x-(R-(R^2-y.^2).^0.5);
% y2=R*asin(y./R);
% figure;
% plot(x2,y2,'g-');
% axis([-1 1 -1 1]);

%%
%¶ÁÈ¡Í¼Æ¬
close all;clc;clear;
% [img_name,img_path]=uigetfile('*.*');
% img_original=imread(fullfile(img_path,img_name));
% img_original=imread('E:\20150412\µØÕð²¨\4.jpg');
% figure;imshow(img_original);hold on;
% if size(img_original,3)>1
%     img_original=double(rgb2gray(img_original));
% end
% initial_points=zeros(2,2);
% for i=1:2
%     initial_points(:,i)=ginput(1)';
%      plot(initial_points(1, i), initial_points(2, i), 'rx');
%     hold on;
% end
first_time=1;
if 1==first_time
    [img_name, img_path]=uigetfile('*.*');
    img=imread(fullfile(img_path, img_name));
    %% get the top-left and bottom-right corner of the rectangle window where we perform TILT.

    figure(1);
    imshow(img);
    hold on;
    initial_points=zeros(2, 2);
    for i=1:2
        initial_points(:, i)=ginput(1)';
        plot(initial_points(1, i), initial_points(2, i), 'rx');
        hold on;
    end
    save data.mat img_name img_path initial_points;
else
    load data.mat img_name img_path initial_points;
     img=imread(fullfile(img_path, img_name));
     figure(1);
    imshow(img);
    hold on;
     for i=1:2
        plot(initial_points(1, i), initial_points(2, i), 'rx');
        hold on;
    end
end
img=img(initial_points(2, 1):initial_points(2, 2),initial_points(1, 1):initial_points(1, 2));
figure,imshow(img);

%%
%²¨ÐÎ½ÃÕý
[h,w]=size(img);
cx=(w+1)/2;
cy=(h+1)/2;
I=zeros(h,w);
I0=img;
R=600;
Svd=[];
Rank1=[];
nu=[];
% for R=1000
for R=300
for x=1:w;
    for y=1:h
        xx=x-cx;
        yy=y-cy;
%         xx=xx+(R-(R^2-yy.^2).^0.5);
%         yy=R*sin(yy./R);
       
        xx=xx-(R-(R^2-yy.^2).^0.5);
        yy=R*asin(yy./R);

        xx = xx + cx;
        yy = yy + cy;
        
        x_1 = floor(xx);
        y_1 = floor(yy);
        x_2 = x_1 + 1;
        y_2 = y_1 + 1;
        x_s = min(max(x_1,1),w);
        y_s = min(max(y_1,1),h);
        f1 = I0(y_s,x_s);
        
        y_s = min(max(y_2,1),h);
        f2 = I0(y_s,x_s);
        
        x_s = min(max(x_2,1),w);
        y_s = min(max(y_1,1),h);
        f3 = I0(y_s,x_s);
        y_s = min(max(y_2,1),h);
        f4 = I0(y_s,x_s);
        f13 = (xx-x_1)*f3+(x_2-xx)*f1;
        f24 = (xx-x_1)*f4+(x_2-xx)*f2;

        I(y, x) = (yy - y_1)*f24+(y_2-yy)*f13;
    end
end
figure(R),imshow(I/255);
Svd=[Svd svd(I)];
Rank1=[Rank1 rank(I)];
nu=[nu sum(svd(I))];
end
Svd=[Svd svd(double(img))];
Rank1=[Rank1 rank(double(img))];
figure,imshow(img/255);
figure,imshow(I/255);
 imwrite(double(I),'curve0430.bmp');















