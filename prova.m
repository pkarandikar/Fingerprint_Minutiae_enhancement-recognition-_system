clc;clear;close all;
cartella_lettura   = 'DB1_B\';

scrivi=1
if scrivi
    for ii=104:110
        for jj=4:8
            close all;
            nome = strcat(cartella_lettura,num2str(ii),'_',num2str(jj),'.tif')
            a = imread(nome);
            a0 = a;
            [oimg,fimg,bwimg,eimg,enhimg] =  fft_enhance_cubs(a);
            a = enhimg;
            [x0,y0]=supercore6(a,a0);
           
            disp('fatto!');
            pause;
           
        end
    end
else
    for ii=102:110
        for jj=1:8
            close all;
            stringa = strcat('pc',num2str(ii),'_',num2str(jj),'.fig')
            hgload(stringa);
            pause;
        end
    end
end
