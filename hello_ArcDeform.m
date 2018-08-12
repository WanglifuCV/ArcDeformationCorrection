clc;
close all;
task_name='TILT_ANGLE';
%% load the image
first_time=1;

%% Debuging Parameters
isChoosePic = 1;
isSimulate = 1;


if 1==first_time
    if isChoosePic == 1
        [img_name, img_path]=uigetfile('*.*');
        img=imread(fullfile(img_path, img_name));
    else
        if isSimulate == 0
            img = imread( '4.jpg' );
            img_name = '4.jpg';
            img_path = 'E:\PhD\地震图\研究和总结\地震波的弧线问题\Code\arc Lifu\';
        elseif isSimulate == 1
            img = imread( 'Wave_Square.tif' );
            img = imresize( img, 0.4 );
            img_name = 'Wave_Square.tif';
            img_path = 'E:\PhD\地震图\研究和总结\地震波的弧线问题\Code\arc Lifu\';
            tau_orig(1) = 300;
            tau_orig(2) = -20;
            tau_orig(2) = tau_orig(2)*pi/180;
            tau_orig(3) = 10;
            img = im2uint8(generate_curved_image( img, tau_orig ));
        elseif isSimulate == 2
            img = imread( 'Curve1.png' );
            img_name = 'Curve1.png';
            img_path = 'E:\PhD\地震图\研究和总结\地震波的弧线问题\Code\arc Lifu\';
        elseif isSimulate == 3
            img = imread( 'Waveform.png' );
%             img = imresize( img, 0.5 );
            img_name = 'Waveform.png';
            img_path = 'E:\PhD\地震图\研究和总结\地震波的弧线问题\Code\arc Lifu\';
            tau_orig(1) = 600;
            tau_orig(2) = -10;
            tau_orig(2) = tau_orig(2)*pi/180;
            tau_orig(3) = 10;
            img = im2uint8(generate_curved_image( img, tau_orig ));
        end
    end
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
        plot(initial_points(1, i), initial_points(2, i), 'rx', 'linewidth', 2 );
        hold on;
    end
end
if size(img,3)>1
    img = double(rgb2gray(img));
end
% if length(size(img))>1
% img=img(:,:,3);
% end
%draw rectangle
width=floor(initial_points(1,2)-initial_points(1,1));
height=floor(initial_points(2,2)-initial_points(2,1));
rectangle('Position',[floor(initial_points(1,1)) floor(initial_points(2,1)) width height], 'linewidth', 2,'edgecolor', 'r' );

%Get 4 Points
%img=double(img);
%[FourPoints]=GetFourPoints(img);
%% Run TILT:
[Dotau, A, E, f,tau, focus_size, error_sign,A_scale]=...
    TILT_ArcDeform(img, 'ANGLE', initial_points, 'SAVE_PATH', fullfile(cd, task_name),...
    'BLUR_SIGMA', 2, 'BLUR_NEIGHBOR', 2,...
    'PYRAMID_MAXLEVEL', 4, 'DISPLAY_INTER', 1);


Dotau = Dotau * ( 255/max(Dotau(:)) );
focus_img = img( round(initial_points(2,1)) : round( initial_points(2,2) ) , round(initial_points(1,1)) : round(initial_points(1,2)) , : );

figure(700);imshow( uint8(Dotau) );
figure(800);imshow( uint8( focus_img ) );

imwrite( uint8(Dotau), 'Rectified.tif' );
imwrite( uint8(focus_img), 'Focus_image.tif' );
save 'Result.mat';