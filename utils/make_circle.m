function Circ = make_circle(mx,my,r)
%function R = make_rect(mx,my,a,b)
%   
%make rectangle function
%input: size of image = [mx,my]
%1/2 length of sides = [a,b]
 
R = zeros(mx,my);
 
Cx = floor(mx/2);
Cy = floor(my/2);
 
[columnsInImage rowsInImage] = meshgrid(1:mx, 1:my);
 
Circ = (rowsInImage - Cy).^2 ...
    + (columnsInImage - Cx).^2 <= r.^2;
 
 
end