function [outx,outy]=supercore(inputimage)


abilitafigura = 0;
abilitaprintf = 0;
fingerprint=inputimage;

img = double(fingerprint);
[dimx,dimy]=size(fingerprint);
mat_ok = sogliavarianza(fingerprint,8,20,10,44);
[oimg] = orientation_image_luigiopt(img);
[gx,gy]=gradient(oimg);
y0=(gx.^2+gy.^2).*double(mat_ok);
y0=wiener2(y0,[30 30]);
y0max = max(max(y0));
y1 = y0>y0max*0.2;
thinned = bwmorph(y1,'thin',Inf);

thinned_ok = zeros(size(thinned));
for ii=1:dimx
    for jj=1:dimy
        if thinned(ii,jj)==1 && ii>1 && ii<dimx && jj>1 && jj<dimy
            blocco = thinned(ii-1:ii+1,jj-1:jj+1);
            if sum(sum(double(blocco)))==2
                thinned_ok(ii,jj)=1;
            end
        end
    end
end
thinned_ok = logical(thinned_ok);
if abilitafigura
    figure('Name','Ending point');
    imshow(thinned_ok);
end
puntidacontrollare = nnz(thinned_ok);
if abilitaprintf
    disp('*********************************');
    disp('Numero punti di ending');
    disp(puntidacontrollare);
end
ending    = thinned_ok;
complesso = filtraggiocomplesso(img,55,1).*double(mat_ok);

bcx = 32;
bcy = 32;
[dimx dimy]=size(ending);
restox = mod(dimx,bcx);
restoy = mod(dimy,bcy);

ending(dimx+restox,dimy+restoy)    = 0;
complesso(dimx+restox,dimy+restoy) = 0;
dimxt = dimx;
dimyt = dimy;
dimx  = dimx+restox;
dimy  = dimy+restoy;

massimoassoluto = max(max(complesso));
soglia  = 0.65;
blocco_massimi = zeros(dimx,dimy);
blocco_logico  = logical(zeros(dimx,dimy));
nbx = dimx/bcx;
nby = dimy/bcy;
for ii=1:nbx
    for jj=1:nby
        x0=(ii-1)*bcx+1;
        x1=ii*bcx;
        y0=(jj-1)*bcy+1;
        y1=jj*bcy;

        blocco = complesso(x0:x1,y0:y1);
        [massimo_vettore,posizione_vettore]=max(blocco);
        [massimo,posizione]=max(massimo_vettore);
        y_max=posizione;
        x_max=posizione_vettore(posizione);


        if massimo>=soglia*massimoassoluto
            posx = x_max+(ii-1)*bcx;
            posy = y_max+(jj-1)*bcy;
            blocco_massimi(posx,posy)        = massimo;
            blocco_logico(posx,posy)         = 1;
        end
    end
end

if abilitafigura
    figure('Name','Blocco > soglia');
    imshow(blocco_logico);
end
if abilitaprintf
    disp('Numero punti dei blocchi > soglia');
    disp(nnz(blocco_logico));
    disp('*********************************');
end
f abilitaprintf
    disp('Per ogni punto del blocco (a partire da massimo)vedo quanti ending ci cadono');
end
[dimx,dimy]=size(ending);
[posx,posy,valore]=find(blocco_massimi);
[valoreordinati,indici]=sort(valore);
posx = posx(indici);
posy = posy(indici);
% finestra = 20;
lunghezza = length(valore);
[distanze,posizione] = bwdist(ending);
if abilitaprintf
    disp(' --- ');
end
distanza_minima  = Inf;
posizione_minima = 0;

for jj=1:lunghezza
    ii    = lunghezza-jj+1;
    if abilitaprintf
        disp('Punto numero...');
        disp(jj);
    end
    cordx = posx(ii);
    cordy = posy(ii);
    if abilitaprintf
        disp('Valore filtro complesso');
        disp(blocco_massimi(cordx,cordy));
    end
    
    distanza_corrente = distanze(cordx,cordy);
    if abilitaprintf
        disp('Distanza da ending piu'' vicino');
        disp(distanza_corrente);
    end
    if distanza_corrente < distanza_minima
        distanza_minima  = distanza_corrente;
        posizione_minima = posizione(cordx,cordy);
        x0 = cordx;
        y0 = cordy;
    end
    if abilitaprintf
        disp(' --- ');
    end
end
x1=0;
y1=0;
min_dis = Inf;
for ii=1:dimx
    for jj=1:dimy
        if ending(ii,jj)
            cor_dis = (ii-x0)^2+(jj-y0)^2;
            if cor_dis < min_dis
                min_dis = cor_dis;
                x1 = ii;
                y1 = jj;
            end
        end
    end
end

outx = x1;
outy = y1;
