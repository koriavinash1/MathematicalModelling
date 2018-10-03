clc
clear all

a=imread('12.png'); % reading the image
a=im2double(a); % normalizing the instensity values to lie between o and 1

ref=a;
noiseL = 80;
ad=imnoise(a,'gaussian', noiseL/255.0); % adding Gaussian noise of mean zero and variance 0.01
timestep=0.2; % timestep size used in numerical approximation
Niter=60; % number of iterations 

alpha=2.7; % Used in Numerical approximation
w= exp(4*alpha/9); % Used in Numerical approximation
b=eed( ad, 0.1, 1, .24, 10, 2, 1); 