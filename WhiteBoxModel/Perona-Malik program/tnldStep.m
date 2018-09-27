function r=tnldStep(L,a,b,c,ip)
% Discrete numerical scheme of dL/dt for tensor diffusion
N=size(L,1);
M=size(L,2);

Lpc=translateImage(L,1,0);
Lpp=translateImage(L,1,1);
Lcp=translateImage(L,0,1);
Lmp=translateImage(L,-1,1);
Lmc=translateImage(L,-1,0);
Lmm=translateImage(L,-1,-1);
Lcm=translateImage(L,0,-1);
Lpm=translateImage(L,1,-1);

amc=translateImage(a,-1,0);
apc=translateImage(a,1,0);
apm=translateImage(a,1,-1);
amp=translateImage(a,-1,1);
app=translateImage(a,1,1);
amm=translateImage(a,-1,-1);
acp=translateImage(a,0,1);
acm=translateImage(a,0,-1);



bmc=translateImage(b,-1,0);
bcm=translateImage(b,0,-1);
bpc=translateImage(b,1,0);
bcp=translateImage(b,0,1);
bmp=translateImage(b,-1,1);
bpp=translateImage(b,1,1);
bmm=translateImage(b,-1,-1);
bpm=translateImage(b,1,-1);

ccp=translateImage(c,0,1);
ccm=translateImage(c,0,-1);
cmc=translateImage(c,-1,0);
cpc=translateImage(c,1,0);
cpm=translateImage(c,1,-1);
cmp=translateImage(c,-1,1);
cpp=translateImage(c,1,1);
cmm=translateImage(c,-1,-1);



% Weickert and Scharr weights
      %J1=a.*(-3/32*Lmp+3/32*Lpp-10/32*Lmc+10/32*Lpc-3/32*Lmm+3/32*Lpm)+b.*(3/32*Lmp+10/32*Lcp+3/32*Lpp-3/32*Lmm-10/32*Lcm-3/32*Lpm);
      %J2=b.*(-3/32*Lmp+3/32*Lpp-10/32*Lmc+10/32*Lpc-3/32*Lmm+3/32*Lpm)+c.*(3/32*Lmp+10/32*Lcp+3/32*Lpp-3/32*Lmm-10/32*Lcm-3/32*Lpm);

% Using new 9 point weights obtained by fpm method
%J1=a.*(-0.0276*Lmp+0.0276*Lpp-0.4447*Lmc+0.4447*Lpc-0.0276*Lmm+0.0276*Lpm)+b.*(0.0276*Lmp+0.4447*Lcp+0.0276*Lpp-0.0276*Lmm-0.4447*Lcm-0.0276*Lpm);
%J2=b.*(-0.0276*Lmp+0.0276*Lpp-0.4447*Lmc+0.4447*Lpc-0.0276*Lmm+0.0276*Lpm)+c.*(0.0276*Lmp+0.4447*Lcp+0.0276*Lpp-0.0276*Lmm-0.4447*Lcm-0.0276*Lpm);
% end of fpm weights

% Sobel filter weights
%J1=a.*(-1/8*Lmp+1/8*Lpp-2/8*Lmc+2/8*Lpc-1/8*Lmm+1/8*Lpm)+b.*(1/8*Lmp+2/8*Lcp+1/8*Lpp-1/8*Lmm-2/8*Lcm-1/8*Lpm);
%J2=b.*(-1/8*Lmp+1/8*Lpp-2/8*Lmc+2/8*Lpc-1/8*Lmm+1/8*Lpm)+c.*(1/8*Lmp+2/8*Lcp+1/8*Lpp-1/8*Lmm-2/8*Lcm-1/8*Lpm);

% Kroon filters
%J1=a.*(-0.006*Lmp+0.006*Lpp-0.948*Lmc+0.948*Lpc-0.006*Lmm+0.006*Lpm)+b.*(0.006*Lmp+0.948*Lcp+0.006*Lpp-0.006*Lmm-0.948*Lcm-0.006*Lpm);
%J2=b.*(-0.006*Lmp+0.006*Lpp-0.948*Lmc+0.948*Lpc-0.006*Lmm+0.006*Lpm)+c.*(0.006*Lmp+0.948*Lcp+0.006*Lpp-0.006*Lmm-0.948*Lcm-0.006*Lpm);
% end of kroon


% J1mc=translateImage(J1,-1,0);
% J1pc=translateImage(J1,1,0);
% J1mp=translateImage(J1,-1,1);
% J1pp=translateImage(J1,1,1);
% J1mm=translateImage(J1,-1,-1);
% J1pm=translateImage(J1,1,-1);
% 
% J2cm=translateImage(J2,0,-1);
% J2cp=translateImage(J2,0,1);
% J2mp=translateImage(J2,-1,1);
% J2pp=translateImage(J2,1,1);
% J2mm=translateImage(J2,-1,-1);
% J2pm=translateImage(J2,1,-1);



% standar scheme
%r=-1/4*(bmc+bcp).*Lmp+1/2*(ccp+c).*Lcp+1/4*(bpc+bcp).*Lpp+1/2*(amc+a).*Lmc-1/2*(amc+2*a+apc+ccm+2*c+ccp).*L+1/2*(apc+a).*Lpc+1/4*(bmc+bcm).*Lmm+1/2*(ccm+c).*Lcm-1/4*(bpc+bcm).*Lpm;

% Non negative scheme

%r=1/4*((abs(bmp)-bmp)+(abs(b)-b)).*Lmp+1/2*(ccp+c-abs(bcp)-abs(b)).*Lcp+1/4*(abs(bpp)+bpp+abs(b)+b).*Lpp+1/2*(amc+a-abs(bmc)-abs(b)).*Lmc+(1/2*(-amc-2*a-apc)-1/4*(abs(bmp)-bmp+abs(bpp)+bpp)-1/4*(abs(bmm)+bmm+abs(bpm)-bpm)+1/2*(abs(bmc)+abs(bpc)+abs(bcm)+abs(bcp)+2*abs(b))-1/2*(ccm+2*c+ccp)).*L+1/2*(apc+a-abs(bpc)-abs(b)).*Lpc+1/4*(abs(bmm)+bmm+abs(b)+b).*Lmm+1/2*(ccm+c-abs(bcm)-abs(b)).*Lcm+1/4*(abs(bpm)-bpm+abs(b)-b).*Lpm;
val=ip;
switch val
    case 1
    % New 5*5 stencil scheme
    % Weickert and Scharr weights
      J1=a.*(-3/32*Lmp+3/32*Lpp-10/32*Lmc+10/32*Lpc-3/32*Lmm+3/32*Lpm)+b.*(3/32*Lmp+10/32*Lcp+3/32*Lpp-3/32*Lmm-10/32*Lcm-3/32*Lpm);
      J2=b.*(-3/32*Lmp+3/32*Lpp-10/32*Lmc+10/32*Lpc-3/32*Lmm+3/32*Lpm)+c.*(3/32*Lmp+10/32*Lcp+3/32*Lpp-3/32*Lmm-10/32*Lcm-3/32*Lpm);
      
      J1mc=translateImage(J1,-1,0);
      J1pc=translateImage(J1,1,0);
      J1mp=translateImage(J1,-1,1);
      J1pp=translateImage(J1,1,1);
      J1mm=translateImage(J1,-1,-1);
      J1pm=translateImage(J1,1,-1);

      J2cm=translateImage(J2,0,-1);
      J2cp=translateImage(J2,0,1);
      J2mp=translateImage(J2,-1,1);
      J2pp=translateImage(J2,1,1);
      J2mm=translateImage(J2,-1,-1);
      J2pm=translateImage(J2,1,-1);
      
      % standar scheme
      %r=-1/4*(bmc+bcp).*Lmp+1/2*(ccp+c).*Lcp+1/4*(bpc+bcp).*Lpp+1/2*(amc+a).*Lmc-1/2*(amc+2*a+apc+ccm+2*c+ccp).*L+1/2*(apc+a).*Lpc+1/4*(bmc+bcm).*Lmm+1/2*(ccm+c).*Lcm-1/4*(bpc+bcm).*Lpm;
      
      % Non negative scheme

      r=1/4*((abs(bmp)-bmp)+(abs(b)-b)).*Lmp+1/2*(ccp+c-abs(bcp)-abs(b)).*Lcp+1/4*(abs(bpp)+bpp+abs(b)+b).*Lpp+1/2*(amc+a-abs(bmc)-abs(b)).*Lmc+(1/2*(-amc-2*a-apc)-1/4*(abs(bmp)-bmp+abs(bpp)+bpp)-1/4*(abs(bmm)+bmm+abs(bpm)-bpm)+1/2*(abs(bmc)+abs(bpc)+abs(bcm)+abs(bcp)+2*abs(b))-1/2*(ccm+2*c+ccp)).*L+1/2*(apc+a-abs(bpc)-abs(b)).*Lpc+1/4*(abs(bmm)+bmm+abs(b)+b).*Lmm+1/2*(ccm+c-abs(bcm)-abs(b)).*Lcm+1/4*(abs(bpm)-bpm+abs(b)-b).*Lpm;
      
      % Weickert 5*5 stencil

      %r=-3/32*J1mp+3/32*J1pp-10/32*J1mc+10/32*J1pc-3/32*J1mm+3/32*J1pm+3/32*J2mp+10/32*J2cp+3/32*J2pp-3/32*J2mm-10/32*J2cm-3/32*J2pm;
    case 2
        % New scheme based on fpm method using first derivative approximation
        
        % Using new 9 point weights obtained by fpm method
        J1=a.*(-0.0276*Lmp+0.0276*Lpp-0.4447*Lmc+0.4447*Lpc-0.0276*Lmm+0.0276*Lpm)+b.*(0.0276*Lmp+0.4447*Lcp+0.0276*Lpp-0.0276*Lmm-0.4447*Lcm-0.0276*Lpm);
        J2=b.*(-0.0276*Lmp+0.0276*Lpp-0.4447*Lmc+0.4447*Lpc-0.0276*Lmm+0.0276*Lpm)+c.*(0.0276*Lmp+0.4447*Lcp+0.0276*Lpp-0.0276*Lmm-0.4447*Lcm-0.0276*Lpm);
        % end of fpm weights
        J1mc=translateImage(J1,-1,0);
        J1pc=translateImage(J1,1,0);
        J1mp=translateImage(J1,-1,1);
        J1pp=translateImage(J1,1,1);
        J1mm=translateImage(J1,-1,-1);
        J1pm=translateImage(J1,1,-1);

        J2cm=translateImage(J2,0,-1);
        J2cp=translateImage(J2,0,1);
        J2mp=translateImage(J2,-1,1);
        J2pp=translateImage(J2,1,1);
        J2mm=translateImage(J2,-1,-1);
        J2pm=translateImage(J2,1,-1);
        
        r=-0.0276*J1mp+0.0276*J1pp-0.4447*J1mc+0.4447*J1pc-0.0276*J1mm+0.0276*J1pm+0.0276*J2mp+0.4447*J2cp+0.0276*J2pp-0.0276*J2mm-0.4447*J2cm-0.0276*J2pm;

        % New 5*5 SOBEL filter scheme

        %r=-1/8*J1mp+1/8*J1pp-2/8*J1mc+2/8*J1pc-1/8*J1mm+1/8*J1pm+1/8*J2mp+2/8*J2cp+1/8*J2pp-1/8*J2mm-2/8*J2cm-1/8*J2pm;

        % Kroon scheme
        %r=-0.006*J1mp+0.006*J1pp-0.948*J1mc+0.948*J1pc-0.006*J1mm+0.006*J1pm+0.006*J2mp+0.948*J2cp+0.006*J2pp-0.006*J2mm-0.948*J2cm-0.006*J2pm;

        % New scheme based on 3*3 stencil using fpm weights
        %r=a.*(0.0553*Lmp-0.1106*Lcp+0.0553*Lpp+0.8894*Lmc-1.7788*L+0.8894*Lpc+0.0553*Lmm-0.1106*Lcm+0.0553*Lpm)+(-0.0276*Lmp+0.0276*Lpp-0.4447*Lmc+0.4447*Lpc-0.0276*Lmm+0.0276*Lpm).*(-0.0276*amp+0.0276*app-0.4447*amc+0.4447*apc-0.0276*amm+0.0276*apm+0.0276*bmp+0.4447*bcp+0.0276*bpp-0.0276*bmm-0.4447*bcm-0.0276*bpm)+(0.0276*Lmp+0.4447*Lcp+0.0276*Lpp-0.0276*Lmm-0.4447*Lcm-0.0276*Lpm).*(-0.0276*bmp+0.0276*bpp-0.4447*bmc+0.4447*bpc-0.0276*bmm+0.0276*bpm+0.0276*cmp+0.4447*ccp+0.0276*cpp-0.0276*cmm-0.4447*ccm-0.0276*cpm)+2*b.*(-0.25*Lmp+0.25*Lpp+0.25*Lmm-0.25*Lpm)+c.*(0.0553*Lmp+0.8894*Lcp+0.0553*Lpp-0.1106*Lmc-1.7788*L-0.1106*Lpc+0.0553*Lmm+0.8894*Lcm+0.0553*Lpm);
        case 3
        % sobel 5*5 stencil scheme
        % Sobel filter weights
        J1=a.*(-1/8*Lmp+1/8*Lpp-2/8*Lmc+2/8*Lpc-1/8*Lmm+1/8*Lpm)+b.*(1/8*Lmp+2/8*Lcp+1/8*Lpp-1/8*Lmm-2/8*Lcm-1/8*Lpm);
        J2=b.*(-1/8*Lmp+1/8*Lpp-2/8*Lmc+2/8*Lpc-1/8*Lmm+1/8*Lpm)+c.*(1/8*Lmp+2/8*Lcp+1/8*Lpp-1/8*Lmm-2/8*Lcm-1/8*Lpm);

      
      J1mc=translateImage(J1,-1,0);
      J1pc=translateImage(J1,1,0);
      J1mp=translateImage(J1,-1,1);
      J1pp=translateImage(J1,1,1);
      J1mm=translateImage(J1,-1,-1);
      J1pm=translateImage(J1,1,-1);

      J2cm=translateImage(J2,0,-1);
      J2cp=translateImage(J2,0,1);
      J2mp=translateImage(J2,-1,1);
      J2pp=translateImage(J2,1,1);
      J2mm=translateImage(J2,-1,-1);
      J2pm=translateImage(J2,1,-1);
      
      % standar scheme
      %r=-1/4*(bmc+bcp).*Lmp+1/2*(ccp+c).*Lcp+1/4*(bpc+bcp).*Lpp+1/2*(amc+a).*Lmc-1/2*(amc+2*a+apc+ccm+2*c+ccp).*L+1/2*(apc+a).*Lpc+1/4*(bmc+bcm).*Lmm+1/2*(ccm+c).*Lcm-1/4*(bpc+bcm).*Lpm;
      
      % Non negative scheme

      %r=1/4*((abs(bmp)-bmp)+(abs(b)-b)).*Lmp+1/2*(ccp+c-abs(bcp)-abs(b)).*Lcp+1/4*(abs(bpp)+bpp+abs(b)+b).*Lpp+1/2*(amc+a-abs(bmc)-abs(b)).*Lmc+(1/2*(-amc-2*a-apc)-1/4*(abs(bmp)-bmp+abs(bpp)+bpp)-1/4*(abs(bmm)+bmm+abs(bpm)-bpm)+1/2*(abs(bmc)+abs(bpc)+abs(bcm)+abs(bcp)+2*abs(b))-1/2*(ccm+2*c+ccp)).*L+1/2*(apc+a-abs(bpc)-abs(b)).*Lpc+1/4*(abs(bmm)+bmm+abs(b)+b).*Lmm+1/2*(ccm+c-abs(bcm)-abs(b)).*Lcm+1/4*(abs(bpm)-bpm+abs(b)-b).*Lpm;
      
      % Weickert 5*5 stencil

      r=-1/8*J1mp+1/8*J1pp-2/8*J1mc+2/8*J1pc-1/8*J1mm+1/8*J1pm+1/8*J2mp+2/8*J2cp+1/8*J2pp-1/8*J2mm-2/8*J2cm-1/8*J2pm;

        otherwise disp('invalid choice')
end


