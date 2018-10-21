function U = eed( L, ref, k, uscale, stepsize, nosteps, verbose, ip, name)
% eed: edge enhancing diffusion

  U = L;
  for i = 1:nosteps
    Ux = gD( U, uscale, 1, 0 );
    Uy = gD( U, uscale, 0, 1 );
    UR2 = Ux.^2 + Uy.^2;
    UR = sqrt(UR2);

    c2 = exp( - (UR / k).^2 );
    c1 = 1/255. * c2;

    d11 = (c1 .* Ux.^2 + c2 .* Uy.^2) ./ (UR2+eps);
    d12 = (c2-c1) .* Ux .* Uy ./ (UR2+eps);
    d22 = (c1 .* Uy.^2 + c2 .* Ux.^2) ./ (UR2+eps);

    U = U + stepsize * tnldStep( U, d11, d12, d22, ip );
  end
  
if verbose
fig = figure(verbose);
subplot(1,2,1); imshow(L,[]); 
title('Original Image');
subplot(1,2,2); imshow(U,[]);
title('Edge Enhancing Diffusion');
% saveas(fig, name);
end
end