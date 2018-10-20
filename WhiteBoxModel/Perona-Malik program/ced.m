function R = ced( L, k, obsscale, intscale, stepsize, nosteps, verbose, ip, name)
% ced: coherence enhancing diffusion
R = L;
for i = 1:nosteps
Rx = gD( R, obsscale, 1, 0 );
Ry = gD( R, obsscale, 0, 1 );
s11 = gD( Rx.^2, intscale, 0, 0 );
s12 = gD( Rx.*Ry, intscale, 0, 0 );
s22 = gD( Ry.^2, intscale, 0, 0 );
alpha = sqrt( (s11-s22).^2 + 4*s12.^2 );
el1 = 1/2 * (s11 + s22 - alpha );
el2 = 1/2 * (s11 + s22 + alpha);

c1 = max(0.01, 1-exp( -(el1-el2).^2 / k^2 ));
c2 = 0;

d11 = 1/2 * (c1+c2+(c2-c1).*(s11-s22)./(alpha+eps));
d12 = (c2-c1).*s12./(alpha+eps);
d22 = 1/2 * (c1+c2-(c2-c1).*(s11-s22)./(alpha+eps));

R = R + stepsize * tnldStep( R, d11, d12, d22, ip );
end

if verbose
fig = figure(verbose);
subplot(1,2,1); imshow(L,[]); 
title('Original Image');
subplot(1,2,2); imshow(R,[]);
title('Coherence Enhancing Diffusion');
saveas(fig, name);
end
end