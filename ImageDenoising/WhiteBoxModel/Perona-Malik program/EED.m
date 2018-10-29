function g=EED(f,ref,k,stepsize,nosteps,verbose,w,ip)
% Perona and Malik diffusion 

 if verbose
    figure(verbose);
    subplot(1,2,1);
    imshow(f);
    title('Original Image');
    drawnow;
 end
 g=f;
 [n,m]=size(f);
 N=n*m;

for i=1:nosteps
    
     g=gD(g,0.001,0,0); % Apply this for Catte et al model
     gx=gD(g,1,1,0);
     gy=gD(g,1,0,1);

   grad2 = gx.*gx + gy.*gy;
   c = EED_D(gx, gy, 'x');
%   Dy = EED_D(gx, gy, 'y');
%   c  = sqrt(Dx .*Dx + Dy .*Dy);
 
   
   
   val=ip;
   switch val
       case 1
           % With  Weickert standard explicit scheme
           g=g+stepsize*snldStep(g,c,w,ip);
          
       otherwise disp('invalid choice');
   end
    
     if verbose
        figure(verbose);
        subplot(1,2,2);
        imshow(g);
        title('Edge Enhanced Diffusion');
        drawnow;
     end
     u=g;    
end