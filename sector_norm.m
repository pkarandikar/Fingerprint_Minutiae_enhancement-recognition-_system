function [disk,vector] = sector_norm( image , mode ,sfasamento)


global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice matricer num_disk


N=h_lato;

M=n_sectors+2;

size_m=N*N;


mean_s=zeros(M,1);
varn_s=zeros(M,1);
num_s=zeros(M,1);
image1=zeros(h_lato);
Mo=100;
Vo=100;

image = double(image);
if sfasamento == 0
    
    if mode==0
        for ( i=1:1:size_m)
            
            tmp=matrice(i);
            tmp=tmp+1;

            mean_s(tmp)= mean_s(tmp)+image(i);
            num_s(tmp)=num_s(tmp)+1;

        end
        for (i=1:1:M)
            mean_s(i)=mean_s(i)/num_s(i);
        end
        for ( i=1:1:size_m)
            tmp=matrice(i);
            
            tmp=tmp+1;

            varn_s(tmp)= varn_s(tmp) + (image(i)- mean_s(tmp))^2;

        end
        for (i=1:1:M)
            varn_s(i)= varn_s(i) / num_s(i);
        end
        for ( i=1:1:size_m)
            tmp=matrice(i);
            
            tmp=tmp+1;
            if (abs(varn_s(tmp))>1)
                if ((image(i) - mean_s(tmp))<0)
                    image1(i)=Mo - (Vo/varn_s(tmp)*((image(i) - mean_s(tmp))^2))^0.5;
                else
                    image1(i)=Mo + (Vo/varn_s(tmp)*((image(i) - mean_s(tmp))^2))^0.5;
                end
            else
                image1(i)=Mo;
            end
        end
        disk=image1;
        vector=varn_s;
    end
   
    if mode==1
        for ( i=1:1:size_m)
            
            tmp=matrice(i);
            tmp=tmp+1;

            mean_s(tmp)= mean_s(tmp)+image(i);
            num_s(tmp)=num_s(tmp)+1;

        end
        for (i=1:1:M)
            mean_s(i)=mean_s(i)/num_s(i);
        end
        for ( i=1:1:size_m)
            tmp=matrice(i);
            
            tmp=tmp+1;

            varn_s(tmp)= varn_s(tmp) + (image(i)- mean_s(tmp))^2;

        end
        for (i=1:1:M)
            varn_s(i)= varn_s(i) / num_s(i);
        end
        for (i=1:1:size_m)
            tmp=matrice(i);
            
            tmp=tmp+1;
            image1(i)=varn_s(tmp);
        end
        vettore=zeros(M,1);
        for ( i=1:1:size_m)
            tmp=matrice(i);
            
            tmp=tmp+1;

            vettore(tmp)= vettore(tmp) + abs(image(i)- mean_s(tmp));

        end
        for (i=1:1:M)
            vettore(i)=vettore(i)/num_s(i);
        end
        disk=image1;
        vector=vettore;
    end
end

if sfasamento == 1
    
    if mode==0
        for ( i=1:1:size_m)
            
            tmp=matricer(i);
            tmp=tmp+1;

            mean_s(tmp)= mean_s(tmp)+image(i);
            num_s(tmp)=num_s(tmp)+1;

        end
        for (i=1:1:M)
            mean_s(i)=mean_s(i)/num_s(i);
        end
        for ( i=1:1:size_m)
            tmp=matricer(i);
            
            tmp=tmp+1;

            varn_s(tmp)= varn_s(tmp) + (image(i)- mean_s(tmp))^2;

        end
        for (i=1:1:M)
            varn_s(i)= varn_s(i) / num_s(i);
        end
        for ( i=1:1:size_m)
            tmp=matricer(i);
            
            tmp=tmp+1;
            if (abs(varn_s(tmp))>1)
                if ((image(i) - mean_s(tmp))<0)
                    image1(i)=Mo - (Vo/varn_s(tmp)*((image(i) - mean_s(tmp))^2))^0.5;
                else
                    image1(i)=Mo + (Vo/varn_s(tmp)*((image(i) - mean_s(tmp))^2))^0.5;
                end
            else
                image1(i)=Mo;
            end
        end
        disk=image1;
        vector=varn_s;
    end
   
    if mode==1
        for ( i=1:1:size_m)
            
            tmp=matricer(i);
            tmp=tmp+1;

            mean_s(tmp)= mean_s(tmp)+image(i);
            num_s(tmp)=num_s(tmp)+1;

        end
        for (i=1:1:M)
            mean_s(i)=mean_s(i)/num_s(i);
        end
        for ( i=1:1:size_m)
            tmp=matricer(i);
            
            tmp=tmp+1;

            varn_s(tmp)= varn_s(tmp) + (image(i)- mean_s(tmp))^2;

        end
        for (i=1:1:M)
            varn_s(i)= varn_s(i) / num_s(i);
        end
        for (i=1:1:size_m)
            tmp=matricer(i);
            
            tmp=tmp+1;
            image1(i)=varn_s(tmp);
        end
        vettore=zeros(M,1);
        for ( i=1:1:size_m)
            tmp=matricer(i);
            
            tmp=tmp+1;

            vettore(tmp)= vettore(tmp) + abs(image(i)- mean_s(tmp));

        end
        for (i=1:1:M)
            vettore(i)=vettore(i)/num_s(i);
        end
        disk=image1;
        vector=vettore;
    end
end