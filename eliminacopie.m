function [out]=eliminacopie(ingresso)


L = size(ingresso,1);
out(1,1)=ingresso(1,1);
out(1,2)=ingresso(1,2);

indice = 2;
if L>=2
    for ii=2:L
        trovato=0;
        for jj=1:ii-1
            if (ingresso(ii,1)==ingresso(jj,1)) && (ingresso(ii,2)==ingresso(jj,2))
                trovato = 1;
                break;
            end
        end
        if trovato==0
            out(indice,1) = ingresso(ii,1);
            out(indice,2) = ingresso(ii,2);
            indice        = indice+1;
        end
    end
end