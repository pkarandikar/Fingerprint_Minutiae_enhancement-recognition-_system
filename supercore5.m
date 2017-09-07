function [outx,outy]=supercore5(inputimage)

soglia_var     = 20;
block_size_var = 8;
se_close       = 10;
se_erosion     = 44;


[dimx,dimy] = size(inputimage);
fingerprint = double(inputimage);

mat_ok            = sogliavarianza(fingerprint,block_size_var,soglia_var,se_close,se_erosion);

[oimg]   = orientation_image_luigiopt(fingerprint);
min_pm   = (oimg<=pi/2);
max_pm   = (oimg>pi/2);
distanze = bwdist(min_pm);
bordo    = ((distanze<=sqrt(2)).*mat_ok & ~min_pm);

complesso_novar0       = filtraggiocomplesso(fingerprint,55,1);
blocco_massimi         = double(bordo).*complesso_novar0;
numero_core_candidati  = nnz(blocco_massimi);
indici_core            = find(blocco_massimi);
valori_core            = blocco_massimi(indici_core);

[val_core_ord,ind_ord] = sort(valori_core(:));
val_core_ord           = flipud(val_core_ord);
ind_ord                = flipud(ind_ord);

[outx0,outy0]          = ind2sub([dimx,dimy],indici_core(ind_ord(1)));

outx = outx0;
outy = outy0;