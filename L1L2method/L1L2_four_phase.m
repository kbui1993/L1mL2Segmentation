%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs grayscale four phase segmentation with L1-alpha*L2 
%TV.
%Input:
%   f: image
%   u1_initial: initialization of u1
%   u2_initial: initialization of u2
%   pm: set of parameters
%
%Output:
%   u1: segmentation result of u1
%   u2: segmentation result of u2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u1,u2,c1,c2,c3,c4] = L1L2_four_phase(f,u1_initial, u2_initial, pm)
    
    %initialize u1 and u2
    u1 = u1_initial;
    u2 = u2_initial;
    
    %compute c1, c2, c3, and c4
    error = 10^(-5)*eye(size(u1));
    c1=sum(f(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
    c2=sum(f(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
    c3=sum(f(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
    c4=sum(f(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
    
    %run segmentation algorithm
    for i = 1:pm.outer_iter
        
        %store old u1 and u2
        u1_old = u1;
        u2_old = u2;
        for j = 1:2
            
            %preinitialize variables for DCA algorithm
            if j == 1
                u = u1;
                r = ((c1-f).^2 - (c2-f).^2-(c3-f).^2+(c4-f).^2).*u2+(c2-f).^2 - (c4-f).^2;
            else
                u = u2;
                r = ((c1-f).^2 - (c2-f).^2-(c3-f).^2+(c4-f).^2).*u1+(c3-f).^2 - (c4-f).^2;
            end
            
            %u update
            u = PDHGLS(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.beta);

            if j == 1
                u1 = u;
            else
                u2 = u;
            end
        end
        
        %update c1, ..., c4
        c1=sum(f(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
        c2=sum(f(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
        c3=sum(f(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
        c4=sum(f(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
        
        %compute relative errors
        relerr1 = norm(u1_old - u1, 'fro')/max([norm(u1_old,'fro'), norm(u1, 'fro'), eps]);
        relerr2 = norm(u2_old - u2, 'fro')/max([norm(u2_old,'fro'), norm(u2, 'fro'), eps]);
        
        %stopping condition
        if relerr1 <1e-4 && relerr2 < 1e-4 && i > 2
            fprintf('Number of DCA inner iterations completed: %d \n \n', i);
            break;
        end
    end
    
    if i == pm.outer_iter
        fprintf('Number of DCA outer iterations completed: %d \n \n', i);
    end
end