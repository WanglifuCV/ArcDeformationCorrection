period=4;
im_width=800;
im_height=400;
img=255*ones(im_height+5,im_width);
for x=1:1:im_width;
y=round((im_height-1)/2*(1+sin(period*2*pi/im_width*x))+1);
img(y:y+5,x)=0;
end
figure,imshow(img/255);
imwrite(img,'sin_curve.bmp');
% x=0:0.1:10*pi;
% figure,plot(x,sin(x),'k');