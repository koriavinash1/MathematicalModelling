clc
clear all

a=imread('tomo.jpg'); % reading the image
a=im2double(a); % normalizing the instensity values to lie between o and 1

ref=a;
ad=imnoise(a,'gaussian', 0.6); % adding Gaussian noise of mean zero and variance 0.01
timestep=0.2; % timestep size used in numerical approximation
Niter=60; % number of iterations 

alpha=2.7; % Used in Numerical approximation
w= exp(4*alpha/9); % Used in Numerical approximation

b=pmc(ad,ref,0.001,timestep,Niter,1,w,1); 
% first argument is the noisy image, 2 is the reference image, 3 is the
% lambda value, 4 is the timestep size, 5 is the no of iterations, 6 is the
% value to show the plot, 7 is the w value used in numerical approximation
% and last argument corresponding to choice of the numerical scheme. 

% b=EED(ad,ref,0.001,timestep,Niter,1,w,1); 