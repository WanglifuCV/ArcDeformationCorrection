%% Load Image
clear;
clc;
close all;

[img_name, img_path]=uigetfile('*.*');
img=imread(fullfile(img_path, img_name));

if img_name ~= 0
    [ Rectified, tau ] = ArcDeformRectification( img );
end

Rectified = floor(uint8(cat( 3, Rectified,Rectified,Rectified )));

imshow( Rectified ),title('Rectified Result'); 