
function radial_filter_bank
close all;
FFTN    =   32;
RSTEPS  =   16;
DELTAR  =   (FFTN/2)/RSTEPS

[x,y]   =   meshgrid(-FFTN/2:FFTN/2-1,-FFTN/2:FFTN/2-1);
r       =   sqrt(x.^2+y.^2);
filter  =   [];
pause;
for r0  =   0:DELTAR:(RSTEPS-1)*DELTAR
     msk = 0.5*(1+cos(pi*(r-r0)/DELTAR));
     msk(abs(r-r0)>DELTAR) =0;
     filter = [filter,msk(:)];
     imagesc(msk);
     pause;
end;
radf    =   filter;
save radial_filters radf; 

fp = fopen('radial_filter.h','w');
for i = 1:size(filter,2)
    k = 1;
    fprintf(fp,'{');
    for j = 1:size(filter,1)
        fprintf(fp,'%f,',filter(j,i));
        if(k == 32) k=0; fprintf(fp,'\n'); end;
        k = k+1;
    end;
    fprintf(fp,'}\n');
end;
fclose(fp);