function R = eed( L, k, uscale, stepsize, nosteps, verbose, ip)
% eed: edge enhancing diffusion

if verbose
figure(verbose);
subplot(1,2,1); imshow(L); title('Original Image'); drawnow;
end
% hRun edge enhancing diusion 'nosteps' times 14b

R = L;
for i = 1:nosteps
Rx = gD( R, uscale, 1, 0 );
Ry = gD( R, uscale, 0, 1 );
Rw2 = Rx.^2 + Ry.^2;
Rw = sqrt(Rw2);

c2 = exp( - (Rw / k).^2 );
c1 = 1/5 * c2


a = (c1 .* Rx.^2 + c2 .* Ry.^2) ./ (Rw2+eps);
b = (c2-c1) .* Rx .* Ry ./ (Rw2+eps);
c = (c1 .* Ry.^2 + c2 .* Rx.^2) ./ (Rw2+eps)

R = R + stepsize * tnldStep( R, a, b, c, ip );

if verbose
figure(verbose);
subplot(1,2,2); imshow(R);
title('Edge Enhancing Diffusion'); drawnow;
end
end