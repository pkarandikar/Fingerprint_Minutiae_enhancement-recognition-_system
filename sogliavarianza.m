function [out] = sogliavarianza(fingerprint,blocksize,soglia,se_close,se_erosion)

bxv           = blocksize;
byv           = blocksize;
soglia_var    = soglia;
dimseclose    = se_close;
dimseerode    = se_erosion;

[dimx,dimy]=size(fingerprint);
img = double(fingerprint);

imgd=double(fingerprint);
dimxr=dimx-mod(dimx,bxv);
dimyr=dimy-mod(dimy,byv);
imgr=imgd(1:dimxr,1:dimyr);
nbx=dimxr/bxv;
nby=dimyr/byv;
mat_var=zeros(dimxr,dimyr);
for ii=1:nbx
    for jj=1:nby
        blocco=imgr((ii-1)*bxv+1:ii*bxv,(jj-1)*byv+1:jj*byv);
        media=sum(sum(blocco))/(bxv*byv);
        varianza=1/(bxv*byv)*sum(sum(abs(media.^2-blocco.^2)));
        mat_var((ii-1)*bxv+1:ii*bxv,(jj-1)*byv+1:jj*byv)=sqrt(varianza);
    end
end
mat_ok=zeros(dimxr,dimyr);
pos=find(mat_var>soglia_var);
mat_ok(pos)=1;
mat_ok(dimx,dimy)=0;


bordo = 20;
mat_temp = zeros(size(mat_ok)+2*bordo);
mat_temp(bordo+1:bordo+dimx,bordo+1:bordo+dimy)=mat_ok;
mat_ok = mat_temp;
mat_ok=imclose(mat_ok,ones(dimseclose));%------------------------- close
mat_ok = mat_ok(bordo+1:bordo+dimx,bordo+1:bordo+dimy);


bordo = 20;
mat_temp = zeros(size(mat_ok)+2*bordo);
mat_temp(bordo+1:bordo+dimx,bordo+1:bordo+dimy)=mat_ok;
mat_ok = mat_temp;
mat_ok=imerode(mat_ok,ones(dimseerode));%------------------------- erode
mat_ok = mat_ok(bordo+1:bordo+dimx,bordo+1:bordo+dimy);



out = mat_ok;