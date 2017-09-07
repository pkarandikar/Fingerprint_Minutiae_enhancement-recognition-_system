
function [oimg] = orientation_image_luigiopt_var(x,Ngiven)

alpha   =   0.3;
N       =   Ngiven;
[dimx dimy] = size(x);

msk     =   fspecial('gaussian',7);
x       =   imfilter(x,msk,'symmetric','same');
x       =   pseudo_matched_filter(x,alpha);

hy      =   -fspecial('sobel');
hx      =   transpose(hy);

gx      =   imfilter(x,hx,'symmetric','same');
gy      =   imfilter(x,hy,'symmetric','same');

gmag    =   sqrt(gx.^2+gy.^2);
theta   =   atan(gy./(gx+1e-5));

numeratore   = (gmag.^2).*sin(2*theta);
denominatore = (gmag.^2).*cos(2*theta);
sopra = zeros(dimx+N,dimy+N);
sotto = zeros(dimx+N,dimy+N);
sopra(N/2+1:dimx+N/2,N/2+1:dimy+N/2) = numeratore;
sotto(N/2+1:dimx+N/2,N/2+1:dimy+N/2) = denominatore;
oimg = zeros(dimx,dimy);

vettore_den    = zeros(dimx,1);
vettore_num    = zeros(dimx,1);
num1           = sum(sum(sopra(1:N+1,1:N+1)));
den1           = sum(sum(sotto(1:N+1,1:N+1)));
vettore_den(1) = den1;
vettore_num(1) = num1;

t   = atan2(num1,(den1+1e-5));
t(t<0) = t(t<0)+2*pi;
t   = 0.5*t; 
oimg(1,1)=t;
for ii=2:dimx
    cordx = ii+N/2;
    vettore_num(ii) = vettore_num(ii-1)-sum(sopra(cordx-N/2,1:N+1))+sum(sopra(cordx+N/2,1:N+1));
    vettore_den(ii) = vettore_den(ii-1)-sum(sotto(cordx-N/2,1:N+1))+sum(sotto(cordx+N/2,1:N+1));
    t               = atan2(vettore_num(ii),vettore_den(ii)+1e-5);
    if t<0
        t = t+2*pi;
    end
    
    t   = 0.5*t; 
    oimg(ii,1)=t;
end

vettore_num_dopo = zeros(size(vettore_num));
vettore_den_dopo = zeros(size(vettore_den));
for jj=2:dimy
    cordy = jj+ N/2;
    for ii=1:dimx
        cordx = ii+ N/2;

        vettore_num_dopo(ii) = vettore_num(ii)-sum(sopra(cordx-N/2:cordx+N/2,cordy-N/2))+sum(sopra(cordx-N/2:cordx+N/2,cordy+N/2));
        vettore_den_dopo(ii) = vettore_den(ii)-sum(sotto(cordx-N/2:cordx+N/2,cordy-N/2))+sum(sotto(cordx-N/2:cordx+N/2,cordy+N/2));
        t                    = atan2(vettore_num_dopo(ii),vettore_den_dopo(ii)+1e-5);
        if t<0
            t = t+2*pi;
        end
     
        t   = 0.5*t; 
        oimg(ii,jj)=t;
    end
    vettore_num = vettore_num_dopo;
    vettore_den = vettore_den_dopo;
end




for i = 1:3
    oimg = smoothen_orientation_image(oimg);
end;