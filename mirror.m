function [out,dimxt,dimyt]=mirror(in,aggiunta)


[dimxt,dimyt]=size(in);

in_memo=in;

extx=aggiunta;
exty=aggiunta;

out=zeros(dimxt+2*extx,dimyt+2*exty);
out(extx+1:extx+dimxt,exty+1:exty+dimyt)=in_memo;

for ii=1:extx
    out(extx-ii+1,:)=out(extx+ii,:);
    out(extx+dimxt+ii,:)=out(extx+dimxt-ii+1,:);
end

for ii=1:exty
    out(:,exty-ii+1)=out(:,exty+ii);
    out(:,exty+dimyt+ii)=out(:,exty+dimyt-ii+1);
end