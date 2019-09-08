function h = rescale_color_image(F)
 
%rescale image values to be between 0 and 1
h = F;
for i = 1:3
    h(:,:,i) = rescale_image(F(:,:,i));
end
end