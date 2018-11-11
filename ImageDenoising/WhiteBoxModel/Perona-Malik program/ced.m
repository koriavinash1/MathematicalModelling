function U = ced( L, ref, k, upscale, scale, stepsize, nosteps, verbose, ip, name)
% ced: coherence enhancing diffusion
ip
name
U = L;
for i = 1:nosteps
Ux = gD( U, upscale, 1, 0 );
Uy = gD( U, upscale, 0, 1 );
Js11 = gD( Ux.^2, scale, 0, 0 );
Js12 = gD( Ux.*Uy, scale, 0, 0 );
Js22 = gD( Uy.^2, scale, 0, 0 );
alpha = sqrt( (Js11-Js22).^2 + 4*Js12.^2 );
el1 = 1/2 * (Js11 + Js22 - alpha );
el2 = 1/2 * (Js11 + Js22 + alpha);

c1 = max(0.01, 1-exp( -(el1-el2).^2 / k^2 ));
c2 = 0;

d11 = 1/2 * (c1+c2+(c2-c1).*(Js11-Js22)./(alpha+eps));
d12 = (c2-c1).*Js12./(alpha+eps);
d22 = 1/2 * (c1+c2-(c2-c1).*(Js11-Js22)./(alpha+eps));

U = U + stepsize * tnldStep( U, d11, d12, d22, ip );
end

if verbose
fig = figure(verbose);
subplot(1,2,1); imshow(L,[]); 
title('Input Image');
subplot(1,2,2); imshow(U,[]);
title('Coherence Enhancing Diffusion');
saveas(fig, name);
end
end