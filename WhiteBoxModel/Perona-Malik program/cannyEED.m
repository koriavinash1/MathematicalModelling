clc
clear all

MI_array = [];
MSE_array = [];
SSIM_array = [];
PSNR_array = [];

images = {"01.png", "02.png", "03.png", "04.png", "05.png", "06.png", "07.png", "08.png", "09.png", "10.png", "11.png", "12.png", "tomo.jpg", "triangle.jpg"};
for i=1:length(images)
  a=imread(images{i}); % reading the image
  a=im2double(a); % normalizing the instensity values to lie between o and 1

  ref=a;
  ad = imnoise(a,'gaussian', 0.); % adding Gaussian noise of mean zero and variance 0.01
  
  timestep=0.2; % timestep size used in numerical approximation
  Niter=60; % number of iterations 

  alpha=2.7; % Used in Numerical approximation
  w= exp(4*alpha/9); % Used in Numerical approximation

  b=ceed(ad,ref,0.001,timestep,Niter,0,w, 1, strcat("CannyEED_", images{i})); 
  MI_array = cat(1, MI_array, MI(ref, b));
  MSE_array = cat(1, MSE_array, MSE(ref, b));
  SSIM_array = cat(1, SSIM_array, ssim_index(ref, b)(1));
  PSNR_array = cat(1, PSNR_array, psnr(ref, b));
end

mi = mean(MI_array)
mse = mean(MSE_array)
ssim = mean(SSIM_array)
psnr_ = mean(PSNR_array)