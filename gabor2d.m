function [out] = gabor2d(ingresso,sigmax,sigmay,omega0,teta,padding)


if nargin == 6
else
    padding = 0;
end

[dimx,dimy] = size(ingresso);
out_temp    = zeros(dimx,dimy);

for jj=1:dimy
    out_temp(:,jj) = gabor1d(ingresso(:,jj),sigmax,omega0*cos(pi+teta),padding);
end
out = zeros(dimx,dimy);
for ii=1:dimx
    out(ii,:)      = gabor1d(out_temp(ii,:),sigmay,omega0*sin(teta),padding);
end

