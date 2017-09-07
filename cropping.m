function [CroppedPrint] = cropping(XofCenter,YofCenter,CentralizedPrint)

global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice matricer num_disk



N = h_lato;
M=size(CentralizedPrint,1);

imgN=size(CentralizedPrint,1);
imgM=size(CentralizedPrint,2);

if (YofCenter-floor(N/2)<1)||(YofCenter+floor(N/2)>imgN)||(XofCenter-floor(N/2)<1)||(XofCenter+floor(N/2)>imgM)
    
    temp=zeros(imgN+2*h_lato,imgM+2*h_lato);
    temp(h_lato+1:h_lato+imgN,h_lato+1:h_lato+imgM)=CentralizedPrint;
    CroppedPrint=temp(YofCenter-floor(N/2)+h_lato:YofCenter+floor(N/2)+h_lato,XofCenter-floor(N/2)+h_lato:XofCenter+floor(N/2)+h_lato);
    return;
else
    CroppedPrint=CentralizedPrint(YofCenter-floor(N/2):YofCenter+floor(N/2),XofCenter-floor(N/2):XofCenter+floor(N/2));
    return;
end