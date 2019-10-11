%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs grayscale two phase segmentation using isotropic
%TV.
%Input:
%   f: image
%   u_initial: initialization
%   pm: set of parameters
%
%Output:
%   u: segmentation result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = isoTV_two_phase(f,u_initial, pm)

    %set u
    u = u_initial;
    
    %compute c1 and c2
    error = 10^(-5)*eye(size(u));
    c1=sum(f(:).*u(:))./sum(u(:)+error(:));
    c2=sum(f(:).*(1-u(:)))./sum((1-u(:))+error(:));
    
    %compute r
    r = (c1-f).^2- (c2-f).^2;
    
    %run DCA algorithm
    for i = 1:pm.outer_iter
        u_old = u;
        u = isoPDHG(u, r, pm.lambda, pm.inner_iter, pm.tau, pm.sigma);
        
        %update c1 and c2
        error = 10^(-5)*eye(size(u));
        c1=sum(f(:).*u(:))./sum(u(:)+error(:));
        c2=sum(f(:).*(1-u(:)))./sum((1-u(:))+error(:));
        
        %update r
        r = (c1-f).^2- (c2-f).^2;
        
        %compute relative error
        relerr = norm(u_old - u)/max([norm(u_old), norm(u), eps]);
        
        %stopping condition
        if relerr <5e-3 && i > 2
            fprintf('Number of iterations completed: %d \n \n', i);
            break;
        end
    end
    
    if i == pm.outer_iter
        fprintf('Number of iterations completed: %d \n \n', i);
    end
end