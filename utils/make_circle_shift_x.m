function Circ = make_circle_shift_x(mx,my,r, shift_x)
%function R = make_rect(mx,my,a,b)
%   
%make rectangle function
%input: size of image = [mx,my]
%1/2 length of sides = [a,b]
 
R = zeros(mx,my);
 
Cx = floor(mx/2)+shift_x;
Cy = floor(my/2);
 
[columnsInImage rowsInImage] = meshgrid(1:mx, 1:my);
 
Circ = (rowsInImage - Cy).^2 ...
    + (columnsInImage - Cx).^2 <= r.^2;
 
 
end