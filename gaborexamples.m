clc;clear;close all;
sigmax = 20;
sigmay = 10;
omega0 = 0.1*pi/2;
teta   = pi/4;
omegax  = omega0*cos(teta);
omegay  = omega0*sin(teta);

lato = 100;
x=[-lato:1:lato];
y=[-lato:1:lato];
ax=1/(sqrt(2*pi)*sigmax)*exp(-x.^2/(2*sigmax^2)).*exp(i*omegax*x);
ay=1/(sqrt(2*pi)*sigmay)*exp(-y.^2/(2*sigmay^2)).*exp(i*omegay*y);

a=zeros(length(x),length(y));
for ii=1:length(x)
    for jj=1:length(y)
        a(ii,jj)=1/(sqrt(2*pi)*sigmax)*exp(-x(ii).^2/(2*sigmax^2)).*exp(i*omegax*x(ii))*...
                   1/(sqrt(2*pi)*sigmay)*exp(-y(jj).^2/(2*sigmay^2)).*exp(i*omegay*y(jj));
    end
end

in = 20*exp(-0.002*(x+40).^2);
filtro = 1/(sqrt(2*pi)*sigmax)*exp(-x.^2/(2*sigmax^2)).*exp(i*omega0*x);%----> 1D Gabor filter
y0 = conv(in',filtro);
px = ((length(filtro)-1)+mod((length(filtro)-1),2))/2;
inizio = 1+px; 
fine   = px+length(in);
y0 = y0(inizio:fine);                     % 1d filtered signal (i.e. the central part of convolution)using a standard technique
y1 = gabor1d(in,sigmax,omega0,20);        % 1d filtered signal using recursive filtering
figure,plot(real(y0)),figure,plot(real(y1));
figure,plot(imag(y0)),figure,plot(imag(y1));

filtro = flipud(ax'*ay);
filtro2= flipud(ax')*flipud(ay);


[X,Y] = meshgrid(-5:.05:5, -5:.05:5);                                
z =  3*exp(-(X-2).^2 - (Y-2).^2)+2*exp(-(X+3).^2 - (Y+1).^2);  %---> input 2D signal    

y0 = conv2(z,filtro,'same');                 % 2d filtered signal (i.e. the central part of convolution)using a standard technique
y1 = gabor2d(z,sigmax,sigmay,omega0,teta,10);% 2d filtered signal using recursive filtering


figure,mesh(real(y0)),figure,mesh(real(y1));

figure,mesh(imag(y0)),figure,mesh(imag(y1));