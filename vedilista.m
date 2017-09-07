clc;clear;close all;
cartella_lettura   = 'DB1_B\';

scrivi=1
if scrivi
    for ii=101:110
        for jj=1:8
            close all;
            nome = strcat(cartella_lettura,num2str(ii),'_',num2str(jj),'.tif')
            a = imread(nome);
            a0 = a;
            [oimg,fimg,bwimg,eimg,enhimg] =  fft_enhance_cubs(a);
            a = enhimg;
            [lista]=supercore7_list(a);
            
            imshow(a0);
            hold on;
            for conta=1:size(lista,1)
                plot(lista(conta,2),lista(conta,1),'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',5)
            end
            hold off;
            
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
