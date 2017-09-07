function [outx,outy]=supercore7(inputimage)

soglia_var     = 20;
block_size_var = 8;
se_close       = 10;
se_erosion     = 44;


[dimxt,dimyt] = size(inputimage);
fingerprint = double(inputimage);

complesso_novar   = filtraggiocomplesso(fingerprint,55,1);

punti_trovati     = 100;
while punti_trovati > 10
   
    mat_ok      = sogliavarianza(fingerprint,block_size_var,soglia_var,se_close,se_erosion);
    complesso   = complesso_novar.*mat_ok;
    
    bcx         = 32;
    bcy         = 32;
    restox      = mod(dimxt,bcx);
    restoy      = mod(dimyt,bcy);

    complesso(dimxt+bcx-restox,dimyt+bcy-restoy)  = 0;
    


    dimx   = dimxt+bcx-restox;
    dimy   = dimyt+bcy-restoy;

    massimoassoluto = max(max(complesso));
    soglia          = 0.65;
    blocco_massimi  = zeros(dimx,dimy);
    blocco_logico   = logical(zeros(dimx,dimy));
    nbx             = dimx/bcx;
    nby             = dimy/bcy;
    for ii=1:nbx
        for jj=1:nby
            x0=(ii-1)*bcx+1;
            x1=ii*bcx;
            y0=(jj-1)*bcy+1;
            y1=jj*bcy;

            blocco                              = complesso(x0:x1,y0:y1);
            [massimo_vettore,posizione_vettore] = max(blocco);
            [massimo,posizione]                 = max(massimo_vettore);
            y_max                               = posizione;
            x_max                               = posizione_vettore(posizione);
            if (massimo >= soglia*massimoassoluto)
                posx                      = x_max+(ii-1)*bcx;
                posy                      = y_max+(jj-1)*bcy;
                blocco_massimi(posx,posy) = massimo;
                blocco_logico(posx,posy)  = 1;
            end
        end
    end
    
    punti_trovati = nnz(blocco_logico);
    
    se_erosion    = se_erosion + 5;
end

[oimg]   = orientation_image_luigiopt_var(fingerprint,16);

sfas0     = -40;
sfas1     = +40;
deltasfas = 10;
indice    = 1;

vettore_core   = [];
for sfasamento = sfas0:deltasfas:sfas1
    
    angolo   = pi/2+deg2rad(sfasamento);
    
    min_pm   = (oimg<=angolo);
    
    distanze = bwdist(min_pm);
    bordo    = ((distanze<=sqrt(2)).*mat_ok & ~min_pm);
    bordo    = bwmorph(bordo,'thin','Inf');
   
    complesso_novar0       = complesso_novar;
    blocco_massimi         = double(bordo).*complesso_novar0;
    numero_core_candidati  = nnz(blocco_massimi);
    indici_core            = find(blocco_massimi);
    valori_core            = blocco_massimi(indici_core);

    [val_core_ord,ind_ord] = sort(valori_core(:));
    val_core_ord           = flipud(val_core_ord);
    ind_ord                = flipud(ind_ord);
    
    if length(ind_ord)>=1
        [outx0,outy0]          = ind2sub([dimxt,dimyt],indici_core(ind_ord(1)));
        
        vettore_core(indice,1) = outx0;
        vettore_core(indice,2) = outy0;
        indice                 = indice+1;
    else
    end
    
end


numero       = size(vettore_core,1);


if numero==0
    m             = complesso_novar.*double(mat_ok);
    [valmax,pmax] = max(m(:));
    [outx,outy]   = ind2sub(size(m),pmax);
    return;
end
liste_elemen = zeros(numero,3);
soglia_dist  = 10;
centri       = 0;
for ii=1:numero
    x1 = vettore_core(ii,1);
    y1 = vettore_core(ii,2);
    
    liste_elemen(ii,1) = x1;
    liste_elemen(ii,2) = y1;
    
    if ii==1
        liste_elemen(ii,3) = 1;
        centri             = centri+1;
    else
        trovato = 0;
        for jj=1:(ii-1)
            x0 = vettore_core(jj,1);
            y0 = vettore_core(jj,2);
            if sqrt((x1-x0)^2+(y1-y0)^2)<=soglia_dist
                trovato  = 1;
                pos_trov = jj;
                break;
            end
        end
        if trovato == 1
            liste_elemen(ii,3) = liste_elemen(pos_trov,3);
        else
            liste_elemen(ii,3) = centri+1;
            centri             = centri+1;
        end
    end
end

conta = zeros(centri,1);
for ii=1:numero
    conta(liste_elemen(ii,3))=conta(liste_elemen(ii,3))+1;
end



centro_medio = zeros(centri,2);
for ii=1:numero
    val = liste_elemen(ii,3);
    centro_medio(val,1)=centro_medio(val,1)+liste_elemen(ii,1);
    centro_medio(val,2)=centro_medio(val,2)+liste_elemen(ii,2);
end
for ii=1:centri
    centro_medio (ii,1)=centro_medio(ii,1)/conta(ii);
    centro_medio (ii,2)=centro_medio(ii,2)/conta(ii);
end
centro_medio = round(centro_medio);


indici = find(conta>=3);
if length(indici)>=1
    xmax   = 0;
    for ii=1:length(indici)
        cur_pos = indici(ii);
        if centro_medio(cur_pos,1)>xmax
            xmax = centro_medio(cur_pos,1);
            vs   = cur_pos;
        end
    end

    xmax = 0;
    for ii=1:numero
        if liste_elemen(ii,3)==vs
            if xmax<liste_elemen(ii,1)
                xmax = liste_elemen(ii,1);
                outx = liste_elemen(ii,1);
                outy = liste_elemen(ii,2);
            end
        end
    end
else
    if max(conta)>=2
        [Y,I]= max(conta);
        outx = centro_medio(I,1);
        outy = centro_medio(I,2);
    else
        
        outx = centro_medio(1,1);
        outy = centro_medio(1,2);
    end
end