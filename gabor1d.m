function [out] = gabor1d(ingresso,sigma,omega0,padding)

if sigma >= 3.556
    q = 0.9804*(sigma-3.556)+2.5091;
end
if (sigma >= 0.5) && (sigma < 3.556)
    q = 0.0561*sigma^2+0.5784*sigma-0.2568;
end
if sigma < 0.5
    disp('Error');
    disp('The Gaussian envelope is too narrow, consequently undersampled and, therefore, to be avoided.');
    out = [];
    return;
end

m0    = 1.16680;
m1    = 1.10783;
m2    = 1.40586;

b0    = 1;
scale = (m0+q)*(m1^2+m2^2+2*m1*q+q^2);
b1    = -q*(2*m0*m1+m1^2+m2^2+(2*m0+4*m1)*q+3*q^2)/scale;
b2    = q^2*(m0+2*m1+3*q)/scale;
b3    = -q^3/scale;

B     = (m0*(m1^2+m2^2)/scale)^2;
fd1   = b1*exp(i*1*omega0);
fd2   = b2*exp(i*2*omega0);
fd3   = b3*exp(i*3*omega0);

bd1   = conj(fd1);
bd2   = conj(fd2);
bd3   = conj(fd3);

modifica                              = size(ingresso,1)>1;
if nargin == 4
    dimx                              = length(ingresso);
    ingressot                         = zeros(1,dimx+2*padding);
    ingressot(padding+1:padding+dimx) = ingresso;
    ingresso                          = ingressot;
end


out      = zeros(1,length(ingresso));
out_temp = zeros(1,length(ingresso));

lunghezza   = length(ingresso);
out_temp(1) = ingresso(1)/(1+fd1+fd2+fd3);
out_temp(2) = ingresso(1)/(1+fd1+fd2+fd3);
out_temp(3) = ingresso(1)/(1+fd1+fd2+fd3);
for ii=4:lunghezza
    out_temp(ii) = ingresso(ii)-(fd1*out_temp(ii-1)+fd2*out_temp(ii-2)+fd3*out_temp(ii-3));
end

out(lunghezza)   = B*out_temp(lunghezza)/(1+bd1+bd2+bd3);
out(lunghezza-1) = B*out_temp(lunghezza)/(1+bd1+bd2+bd3);
out(lunghezza-2) = B*out_temp(lunghezza)/(1+bd1+bd2+bd3);
for ii=lunghezza-3:-1:1
    out(ii) = B*out_temp(ii)-(bd1*out(ii+1)+bd2*out(ii+2)+bd3*out(ii+3));
end

if modifica
    out = out';
end
if nargin == 4
    out = out(padding+1:dimx+padding);
end



