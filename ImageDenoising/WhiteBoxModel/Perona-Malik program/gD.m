function g=gD(f,scale,ox,oy)
% Gausssian (Derivative ) convolution
% what's the significance of 3*sigma for kernel size ...
K=ceil(3*scale);
x=-K:K;
Gs=exp(-x.^2/(2*scale^2));
Gs=Gs/sum(Gs);
Gsx=gDerivative(ox,x,Gs,scale);
Gsy=gDerivative(oy,x,Gs,scale);
g=convSepBrd(f,Gsx,Gsy);