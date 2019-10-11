%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function adds gaussian noise to an image.
%Input:
%   I: image to add noise to
%   sigma: standard deviation of gaussian noise
%Output:
%   In: noisy image result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function In = add_noise2(I,sigma)
 
%get image size
[m,n]=size(I);

%convert to double
I = double(I);
 
%Adjust intensities in image I to range from 0 to 1
I = I - min(I(:));
I = I / max(I(:));

%set rng
rng('default');
rng(2);

%create noisy matrix
eta =sigma*randn(m,n);

%add noise
In = I+eta;
 
end