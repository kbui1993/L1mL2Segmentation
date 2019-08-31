%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs the discrete, transpose gradient operator in 
%y-direction.
%
%Input:
%   u: function whose Dy^t will be computed
%Output:
%   d: Dy^t of u
%Note: Code is from Yifei Lou.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = Dyt(u)

%obtain size of u
[rows,cols] = size(u);

%preinitialize
d = zeros(rows,cols);

%compute Dy^t(u)
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
return
