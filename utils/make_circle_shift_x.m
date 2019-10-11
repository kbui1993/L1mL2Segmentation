%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function creates an circle as an initialization for segmentation.
%
%Input:
%   mx: height of the image
%   my: length of the image
%   r: radius of the circle
%   shift_x: shift center of the circle (corresponding to center of image)
%   by shift_x
%Output:
%   Circ: circle region initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Circ = make_circle_shift_x(mx,my,r, shift_x)

%determine center of circle
Cx = floor(mx/2)+shift_x;
Cy = floor(my/2);

%create circle region
[columnsInImage rowsInImage] = meshgrid(1:mx, 1:my);
 
Circ = (rowsInImage - Cy).^2 ...
    + (columnsInImage - Cx).^2 <= r.^2;
 
end