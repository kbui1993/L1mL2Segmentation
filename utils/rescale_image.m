%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function rescales an image so that the image values are between 0 and
%1.
%Input:
%   F: image to be rescaled
%Output:
%   h: rescaled image F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = rescale_image(F)
 
%rescale image values to be between 0 and 1
h = double(F);
h= h-min(min(h));
h = h/max(max(h));
end