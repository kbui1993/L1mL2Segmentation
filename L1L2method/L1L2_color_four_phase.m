%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs color four phase segmentation using L1-alpha*L2 TV.
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
function [u1,u2, c1, c2, c3, c4] = L1L2_color_four_phase(f,u1_initial, u2_initial, pm)

    %separate the channels of f
    f_r = f(:,:,1);
    f_g = f(:,:,2);
    f_b = f(:,:,3);
    
    %initialize u1 and u2
    u1 = u1_initial;
    u2 = u2_initial;
    
    %compute c1, c2, c3, and c4
    error = 10^(-5)*eye(size(u1));
    c1=[sum(f_r(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
        sum(f_g(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
        sum(f_b(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:))];
    c2=[sum(f_r(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
        sum(f_g(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
        sum(f_b(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:))];
    c3=[sum(f_r(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
        sum(f_g(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
        sum(f_b(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:))];
    c4=[sum(f_r(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
        sum(f_g(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
        sum(f_b(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:))];
    
    %run segmentation algorithm
    for i = 1:pm.outer_iter
        
        %store old u1 and u2
        u1_old = u1;
        u2_old = u2;
        for j = 1:2
            
            %preinitialize variables for DCA algorithm
            if j == 1
                u = u1;
                r_r = ((c1(1)-f_r).^2 - (c2(1)-f_r).^2-(c3(1)-f_r).^2+(c4(1)-f_r).^2).*u2+(c2(1)-f_r).^2 - (c4(1)-f_r).^2;
                r_g = ((c1(2)-f_g).^2 - (c2(2)-f_g).^2-(c3(2)-f_g).^2+(c4(2)-f_g).^2).*u2+(c2(2)-f_g).^2 - (c4(2)-f_g).^2;
                r_b = ((c1(3)-f_b).^2 - (c2(3)-f_b).^2-(c3(3)-f_b).^2+(c4(3)-f_b).^2).*u2+(c2(3)-f_b).^2 - (c4(3)-f_b).^2;
                r = r_r+r_g+r_b;
            else
                u = u2;
                r_r = ((c1(1)-f_r).^2 - (c2(1)-f_r).^2-(c3(1)-f_r).^2+(c4(1)-f_r).^2).*u1+(c3(1)-f_r).^2 - (c4(1)-f_r).^2;
                r_g = ((c1(2)-f_g).^2 - (c2(2)-f_g).^2-(c3(2)-f_g).^2+(c4(2)-f_g).^2).*u1+(c3(2)-f_g).^2 - (c4(2)-f_g).^2;
                r_b = ((c1(3)-f_b).^2 - (c2(3)-f_b).^2-(c3(3)-f_b).^2+(c4(3)-f_b).^2).*u1+(c3(3)-f_b).^2 - (c4(3)-f_b).^2;
                r = r_r+r_g+r_b;
            end
            
            %u update
            if strcmp(pm.method, 'PDHG')
                u = PDHG(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.sigma);
            elseif strcmp(pm.method, 'aPDHG')
                u = aPDHG(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.sigma);
            else
                u = Split_Bregman(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.beta);
            end
            if j == 1
                u1 = u;
            else
                u2 = u;
            end
        end
        
        %update c1, ..., c4
        c1=[sum(f_r(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
        sum(f_g(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
        sum(f_b(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:))];
        c2=[sum(f_r(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
        sum(f_g(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
        sum(f_b(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:))];
        c3=[sum(f_r(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
        sum(f_g(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
        sum(f_b(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:))];
        c4=[sum(f_r(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
        sum(f_g(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
        sum(f_b(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:))];
        
        %compute relative errors
        relerr1 = norm(u1_old - u1)/max([norm(u1_old), norm(u1), eps]);
        relerr2 = norm(u2_old - u2)/max([norm(u2_old), norm(u2), eps]);
        
        %stopping condition
        if relerr1 <5e-3 && relerr2 < 5e-3 && i > 2
            fprintf('Number of iterations completed: %d \n \n', i);
            break;
        end
    end
    
    if i == pm.outer_iter
        fprintf('Number of iterations completed: %d \n \n', i);
    end
end