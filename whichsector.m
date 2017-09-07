function [sector_num] = whichsector(index,sfasamento)


global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice matricer num_disk

if sfasamento == 1
    faseaggiunta = deg2rad(180/(num_disk*2));
else
    faseaggiunta = 0;
end

length = h_lato;
x = rem( index , length );
y = floor(index / length);

x = x - floor(length / 2);
y = y - floor(length / 2);

rad = (x*x) + (y*y);
if rad < (h_radius*h_radius)                  
    sector_num = (n_sectors-1)+1;
    sector_num;
    return
end

if rad >= (h_bands*n_bands+h_radius)^2        
    sector_num = (n_sectors-1)+2;
    sector_num;
    return
end   

if x ~= 0
    theta = atan( y / x ) + faseaggiunta;
else 
    if y > 0
        theta = pi/2 + faseaggiunta;
    else
        theta = -pi/2 + faseaggiunta;
    end
end   

if x < 0
    theta = theta + pi;
else
    if theta < 0
        theta = theta + 2*pi;
    end
end

if theta < 0
    theta = theta + 2*pi;
end

r = floor(rad ^ 0.5);
ring = floor(( r-h_radius )/h_bands);
arc = floor(theta /(2*pi/n_arcs));

sector_num = ring * n_arcs + arc;