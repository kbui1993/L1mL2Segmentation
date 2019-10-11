%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function creates an circle as an initialization for segmentation.
%
%Input:
%   mx: height of the image
%   my: length of the image
%   r: radius of the circle
%Output:
%   Circ: circle region initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Circ = make_circle(mx,my,r)

%determine center of circle 
Cx = floor(mx/2);
Cy = floor(my/2);

%create circle
[columnsInImage rowsInImage] = meshgrid(1:mx, 1:my);
 
Circ = (rowsInImage - Cy).^2 ...
    + (columnsInImage - Cx).^2 <= r.^2;
 
 
end