
function angular_filter_bank(BW,fname)
close all;

FFTN    =   32;
TSTEPS  =   12; 
DELTAT  =   pi/TSTEPS;

[x,y]   =   meshgrid(-FFTN/2:FFTN/2-1,-FFTN/2:FFTN/2-1);
r       =   sqrt(x.^2+y.^2);
th      =   atan2(y,x);
th(th<0)=   th(th<0)+2*pi;  
filter  =   [];


for t0  =   0:DELTAT:(TSTEPS-1)*DELTAT
     t1     = t0+pi;   
     d          = angular_distance(th,t0);
     msk        = 1+cos(d*pi/BW); 
     msk(d>BW)  = 0;
     rmsk       = msk;   
     d          = angular_distance(th,t1);
     msk        = 1+cos(d*pi/BW); 
     msk(d>BW)  = 0;
     rmsk       = (rmsk+msk);

     imagesc(rmsk);pause;
     rmsk   = transpose(rmsk);
     filter = [filter,rmsk(:)];
end;
     
     angf  = filter;
     eval(sprintf('save %s angf',fname));


fp = fopen(sprintf('%s.h',fname),'w');
fprintf(fp,'{\n');
for i = 1:size(filter,2)
    i
    k = 1;
    fprintf(fp,'{');
    for j = 1:size(filter,1)
        fprintf(fp,'%f,',filter(j,i));
        if(k == 32) k=0; fprintf(fp,'\n'); end;
        k = k+1;
    end;
    fprintf(fp,'},\n');
end;
fprintf(fp,'};\n');
fclose(fp);

function d = angular_distance(th,t0)
    d = abs(th-t0);
    d = min(d,2*pi-d);