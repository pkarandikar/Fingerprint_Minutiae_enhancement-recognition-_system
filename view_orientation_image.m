
function view_orientation_image(o)
    [h,w]   =   size(o);
    x       =   0:w-1;
    y       =   0:h-1;
    quiver(x,y,cos(o),sin(o)); 
    axis([0 w 0 h]),axis image, axis ij;