%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs color two phase segmentation
%Input:
%   f: image
%   u_initial: initialization
%   pm: set of parameters
%
%Output:
%   u: segmentation result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = L1L2_color_two_phase(f,u_initial, pm)
    %separate the channels of f
    f_r = f(:,:,1);
    f_g = f(:,:,2);
    f_b = f(:,:,3);
    
    %preinitialize u
    u = u_initial;
    
    %compute c1 and c2
    error = 10^(-5)*eye(size(u));
    c1=[sum(f_r(:).*u(:))./sum(u(:)+error(:)); sum(f_g(:).*u(:))./sum(u(:)+error(:)); sum(f_b(:).*u(:))./sum(u(:)+error(:))];
    c2=[sum(f_r(:).*(1-u(:)))./sum((1-u(:))+error(:)); sum(f_g(:).*(1-u(:)))./sum((1-u(:))+error(:)); sum(f_b(:).*(1-u(:)))./sum((1-u(:))+error(:))];
    
    %compute r
    r_r = (c1(1)-f_r).^2- (c2(1)-f_r).^2;
    r_g = (c1(2)-f_g).^2 - (c2(2)-f_g).^2;
    r_b = (c1(3)-f_b).^2 - (c2(3)-f_b).^2;
    r = r_r+r_g+r_b;
    
    %set u
    u = u_initial;
    
    %run segmentation algorithm
    for i = 1:pm.outer_iter
        
        %store old u
        u_old = u;
        
        %u update
        if strcmp(pm.method, 'PDHG')
            u = PDHG(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.sigma);
        elseif strcmp(pm.method, 'aPDHG')
            u = aPDHG(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.sigma);
        else
            u = Split_Bregman(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.beta);
        end
        
        %update c
        c1=[sum(f_r(:).*u(:))./sum(u(:)+error(:)); sum(f_g(:).*u(:))./sum(u(:)+error(:)); sum(f_b(:).*u(:))./sum(u(:)+error(:))];
        c2=[sum(f_r(:).*(1-u(:)))./sum((1-u(:))+error(:)); sum(f_g(:).*(1-u(:)))./sum((1-u(:))+error(:)); sum(f_b(:).*(1-u(:)))./sum((1-u(:))+error(:))];
        
        %update r
        r_r = (c1(1)-f_r).^2- (c2(1)-f_r).^2;
        r_g = (c1(2)-f_g).^2 - (c2(2)-f_g).^2;
        r_b = (c1(3)-f_b).^2 - (c2(3)-f_b).^2;
        r = r_r+r_g+r_b;
        
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