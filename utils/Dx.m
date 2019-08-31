%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs the discrete gradient operator in x-direction.
%
%Input:
%   u: function whose Dx will be computed
%Output:
%   d: Dx of u
%Note: Code is from Yifei Lou.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = Dx(u)
%obtain matrix size of u
[rows,cols] = size(u);

%preinitialize matrix
d = zeros(rows,cols);

%compute Dx(u)
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return