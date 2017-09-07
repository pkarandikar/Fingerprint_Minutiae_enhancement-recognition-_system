function [out] = filtraggiocomplesso(fingerprint,cf_varianza,cf_ordine)


x=[-200:1:200];
y=[-200:1:200];
dimx=size(x,2);
dimy=size(y,2);

varianza   = sqrt(cf_varianza);
ordine     = cf_ordine;
gamma=2;
filtro_core=zeros(dimx,dimy);

for ii=1:dimx
    for jj=1:dimy
        esponente=exp(-(x(ii)^2+y(jj)^2)/(2*varianza^2));
        
        fattore=x(ii)+i*y(jj);
        filtro_core(ii,jj)=esponente*fattore^ordine;
       
    end
end

img=fingerprint;
img=double(img);

[gx,gy]=gradient(img);
num=(gx+i*gy).^2;
den=abs((gx+i*gy).^2);
pos=find(den);
num(pos)=num(pos)./den(pos);
z=zeros(size(img,1),size(img,2));
z(pos)=num(pos);
pos=find(den==0);
z(pos)=1;

temp=z;

[temp,dimxt,dimyt]=mirror(temp,20);
z_f=conv2fft(temp,filtro_core,'same');

z_f=recrop(z_f,dimxt,dimyt,20);
z_f=abs(z_f);

out = z_f;


