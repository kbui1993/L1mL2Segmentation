%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs soft-thresholding
%Input:
%   x: vector/matrix to be thresholded
%   lambda: thresholding parameter
%Output:
%   result: thresholded x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = soft_threshold(x, lambda)
    %performs soft thresholding
    result = sign(x).*max(abs(x) - lambda,0);
end