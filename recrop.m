function [out]=recrop(in,dimxt,dimyt,aggiunta)


extx=aggiunta;
exty=aggiunta;

out=in(extx+1:extx+dimxt,exty+1:exty+dimyt);