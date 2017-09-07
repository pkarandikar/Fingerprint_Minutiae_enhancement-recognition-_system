
a=imread('19_7.bmp');

img=double(a);

[oimg,fimg,bwimg,eimg,enhimg] =  fft_enhance_cubs(img);
img = enhimg;

[xc,yc]=supercore7(img); 

figure('Name','Input fingerprint and core point');
imshow(a);
hold on;
plot(yc,xc,'O');
hold off;