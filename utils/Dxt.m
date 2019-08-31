%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs the discrete, transpose gradient operator in 
%x-direction.
%
%Input:
%   u: function whose Dx^t will be computed
%Output:
%   d: Dx^t of u
%Note: Code is from Yifei Lou.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = Dxt(u)

%obtain size of u
[rows,cols] = size(u);

%preinitialize
d = zeros(rows,cols);

%compute Dx^t(u)
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
return