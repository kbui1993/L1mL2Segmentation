%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs color two phase segmentation with L1-alpha*L2 TV
%Input:
%   f: image
%   u_initial: initialization
%   pm: set of parameters
%
%Output:
%   u: segmentation result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, c1, c2] = L1L2_color_two_phase(f,u_initial, pm)
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
    
    %run segmentation algorithm
    for i = 1:pm.outer_iter
        
        %store old u
        u_old = u;
        
        %u update
        u = PDHGLS(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.beta);
        
        %update c
        c1=[sum(f_r(:).*u(:))./sum(u(:)+error(:)); sum(f_g(:).*u(:))./sum(u(:)+error(:)); sum(f_b(:).*u(:))./sum(u(:)+error(:))];
        c2=[sum(f_r(:).*(1-u(:)))./sum((1-u(:))+error(:)); sum(f_g(:).*(1-u(:)))./sum((1-u(:))+error(:)); sum(f_b(:).*(1-u(:)))./sum((1-u(:))+error(:))];
        
        %update r
        r_r = (c1(1)-f_r).^2- (c2(1)-f_r).^2;
        r_g = (c1(2)-f_g).^2 - (c2(2)-f_g).^2;
        r_b = (c1(3)-f_b).^2 - (c2(3)-f_b).^2;
        r = r_r+r_g+r_b;
        
        %compute relative error
        relerr = norm(u_old - u, 'fro')/max([norm(u_old, 'fro'), norm(u, 'fro'), eps]);
        
        %stopping condition
        if relerr <1e-6 && i > 2
            fprintf('Number of DCA inner iterations completed: %d \n \n', i);
            break;
        end
    end
    
    if i == pm.outer_iter
        fprintf('Number of DCA outer iterations completed: %d \n \n', i);
    end
end