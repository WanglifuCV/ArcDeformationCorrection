function [ Rectified, tau ] = ArcDeformRectification( img )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

task_name='TILT_ANGLE';

figure(1);
imshow(img);
hold on;
initial_points=zeros(2, 2);
for i=1:2
    initial_points(:, i)=ginput(1)';
    plot(initial_points(1, i), initial_points(2, i), 'rx');
     hold on;
end

if size( img, 3 ) > 1
    img = rgb2gray( img );
end

figure(100);imshow( uint8( img( initial_points( 2, 1 ) : initial_points( 2, 2 ), initial_points( 1, 1 ) : initial_points( 1, 2 ) ) ) );
title( 'Input' );

[Dotau, ~, ~, ~,tau, ~, ~,~]=...
    TILT_ArcDeform(img, 'ANGLE', initial_points, 'SAVE_PATH', fullfile(cd, task_name),...
    'BLUR_SIGMA', 2, 'BLUR_NEIGHBOR', 2,...
    'PYRAMID_MAXLEVEL', 4, 'DISPLAY_INTER', 1);

%Dotau = Dotau * ( 255/max(Dotau(:)) );
img=img(initial_points( 2, 1 ) : initial_points( 2, 2 ),1:size(img,2));
[ H, W ] = size( img );

h0 = ( H - 1 )/2;
[ X0, Y0 ] = meshgrid( 1 : W, -h0 : h0 );
XII = X0 + tau (1)*( 1 - cos( (Y0) / tau(1) ) ) - tau(1)*sin( (Y0)/tau(1) )*sin(tau(2));
YII =tau(1)*sin( (Y0)/tau(1) )*cos(tau(2))+tau(3);
img = interp2( X0, Y0, double(img), XII, YII, 'cubic' );
img ( isnan(img ) ) = 255;
Dotau = uint8(img * ( 255/max(img(:)) ));


Rectified = Dotau;


end

