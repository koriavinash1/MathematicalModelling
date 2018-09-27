function r=snldStep(L,c,w,ip)
% Discrete numerical scheme of dL/dt for scalar diffusion
N=size(L,1);
M=size(L,2);
cpc=translateImage(c,1,0);
cmc=translateImage(c,-1,0);
ccp=translateImage(c,0,1);
ccm=translateImage(c,0,-1);


Lpc=translateImage(L,1,0);
Lppc=translateImage(L,2,0);
Lpp=translateImage(L,1,1);
Lcp=translateImage(L,0,1);
Lcpp=translateImage(L,0,2);
Lmp=translateImage(L,-1,1);
Lmc=translateImage(L,-1,0);
Lmmc=translateImage(L,-2,0);
Lmm=translateImage(L,-1,-1);
Lcm=translateImage(L,0,-1);
Lcmm=translateImage(L,0,-2);
%Lcmm=translateImage(L,0,-2);
Lpm=translateImage(L,1,-1);


Lx=-1/(2*(w+2))*Lmp+1/(2*(w+2))*Lpp-w/(2*(w+2))*Lmc+w/(2*(w+2))*Lpc-1/(2*(w+2))*Lmm+1/(2*(w+2))*Lpm;
Ly=1/(2*(w+2))*Lmp+w/(2*(w+2))*Lcp+1/(2*(w+2))*Lpp-1/(2*(w+2))*Lmm-w/(2*(w+2))*Lcm-1/(2*(w+2))*Lpm;


Cpc=C(((Lcp-Lpm+Lcp-Lcm)/4).^2+(Lpc-L).^2);
Cmc=C(((Lcp-Lcm+Lmp-Lmm)/4).^2+(L-Lmc).^2);
Ccp=C(((Lpp-Lmp+Lpc-Lmc)/4).^2+(Lcp-L).^2);
Ccm=C(((Lpc-Lmc+Lpc-L)/4).^2+(L-Lcm).^2);

alpha=0.5;
Cmm=C((alpha/2)*((L-Lmc).^2+(Lcm-Lmm).^2+(L-Lcm).^2+(Lmc-Lmm).^2)+((1-alpha)/2)*((L-Lmm).^2+(Lcm-Lmc).^2));
Cmp=C((alpha/2)*((Lcp-Lmp).^2+(L-Lmc).^2+(Lcp-L).^2+(Lmp-Lmc).^2)+((1-alpha)/2)*((Lmp-L).^2+(Lcp-Lmc).^2));
Cpm=C((alpha/2)*((Lpc-L).^2+(Lpm-Lcm).^2+(Lpc-Lpm).^2+(L-Lcm).^2)+((1-alpha)/2)*((Lpc-Lcm).^2+(Lpm-L).^2));
Cpp=C((alpha/2)*((Lpp-Lcp).^2+(Lpc-L).^2+(Lpp-Lpc).^2+(Lcp-L).^2)+((1-alpha)/2)*((Lpp-L).^2+(Lpc-Lcp).^2));




betap=Ccp+Ccm+Cpc+Cmc;
betad=Cpp+Cmm+Cpm+Cmp;


% For D.Chen scheme put
 lambp=0.5;
 lambd=0.5;
 
val=ip;
switch val
% case 1: weickert  standard scheme or weickert 5*5 stencil scheme
% case 2: new 3*3 scheme or 5*5 scheme based on fpm
% case 3: sobel scheme
    case 1
        r=((cpc+c).*(Lpc-L)-(c+cmc).*(L-Lmc)+(ccp+c).*(Lcp-L)-(c+ccm).*(L-Lcm))/2;
        
              
        
% Scheme acc'ing to aubert
    case 2
        r=lambp.*(Cpc.*Lpc+Cmc.*Lmc+Ccp.*Lcp+Ccm.*Lcm-betap.*L)+0.5*lambd.*(Cpp.*Lpp+Cmm.*Lmm+Cmp.*Lmp+Cpm.*Lpm-betad.*L);
        
    otherwise disp('invalid choice')
end
