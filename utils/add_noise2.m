function In = add_noise2(I,sigma)
 
% function adds gaussian white noise to an image of presribed mean
% and standard deviation sigma
% Fred Park 
% Whittier College
% 9/2014
 
%rand(2);
%rng('default');
%rng(1);
 
[m,n]=size(I);
 
I = double(I);
 
%// Adjust intensities in image I to range from 0 to 1
I = I - min(I(:));
I = I / max(I(:));
 
%// Add noise to image
%v = var(I(:)) / 10^(SNR/10);
%sigma = sqrt(v);
mean_f = 0;
 
rng('default');
rng(2);
eta = mean_f + sigma*randn(m,n);
 
In = I+eta;
 
%In = In - min(I(:));
%In = In / max(I(:));
%In = 255*In;
 
%eta = mean + sigma*randn(m,n); %create noise via randn
%fn = f + eta;