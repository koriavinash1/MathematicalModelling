function g=ceed(f,ref,k,stepsize,nosteps,verbose,w,ip, name)
% Perona and Malik diffusion 
 g=f + edge(ref, 'Canny');
 [n,m]=size(f);
 N=n*m;

for i=1:nosteps
    
     g=gD(g,0.001,0,0); % Apply this for Catte et al model
     gx=gD(g,1,1,0);
     gy=gD(g,1,0,1);

    gpc=translateImage(g,1,0);
    gmc=translateImage(g,-1,0);
    gcp=translateImage(g,0,1);
    gcm=translateImage(g,0,-1);
    gpp=translateImage(g,1,1);
    gmp=translateImage(g,-1,1);
    gpm=translateImage(g,1,-1);
    gmm=translateImage(g,-1,-1);
    size(gmm);
   
   grad2=gx.*gx+gy.*gy;
   c=C(grad2);  
   val=ip;
   switch val
       case 1
           % With  Weickert standard explicit scheme
           g=g+stepsize*snldStep(g,c,w,ip);
          
       otherwise disp('invalid choice');
   end
   u=g;    
end

if verbose
fig = figure(verbose);
subplot(1,2,1); imshow(f,[]); 
title('Input Image');
subplot(1,2,2); imshow(g,[]);
title('Canny EED');
saveas(fig, name);
end
end
