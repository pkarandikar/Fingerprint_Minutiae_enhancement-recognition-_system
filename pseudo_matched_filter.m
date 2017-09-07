
function y = pseudo_matched_filter(x,alpha)
    [h,w]   =   size(x);
    x       =   double(x);
    f       =   fft2(x);
    f       =   (abs(f).^alpha).*f; 
    y       =   real(ifft2(f));