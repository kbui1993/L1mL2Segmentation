function [u1,u2] = L1L2_four_phase(f,u1_initial, u2_initial, pm)
    
    %initialize
    u1 = u1_initial;
    u2 = u2_initial;
    error = 10^(-5)*eye(size(u1));
    c1=sum(f(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
    c2=sum(f(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
    c3=sum(f(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
    c4=sum(f(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
    
    %run DCA algorithm
    for i = 1:pm.outer_iter
        u1_old = u1;
        u2_old = u2;
        for j = 1:2
            if j == 1
                u = u1;
                r = ((c1-f).^2 - (c2-f).^2-(c3-f).^2+(c4-f).^2).*u2+(c2-f).^2 - (c4-f).^2;
            else
                u = u2;
                r = ((c1-f).^2 - (c2-f).^2-(c3-f).^2+(c4-f).^2).*u1+(c3-f).^2 - (c4-f).^2;
            end
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
        c1=sum(f(:).*u1(:).*u2(:))./sum(u1(:).*u2(:)+error(:));
        c2=sum(f(:).*u1(:).*(1-u2(:)))./sum(u1(:).*(1-u2(:))+error(:));
        c3=sum(f(:).*(1-u1(:)).*u2(:))./sum((1-u1(:)).*u2(:)+error(:));
        c4=sum(f(:).*(1-u1(:)).*(1-u2(:)))./sum((1-u1(:)).*(1-u2(:))+error(:));
        relerr1 = norm(u1_old - u1)/max([norm(u1_old), norm(u1), eps]);
        relerr2 = norm(u2_old - u2)/max([norm(u2_old), norm(u2), eps]);
        if relerr1 <5e-3 && relerr2 < 5e-3 && i > 2
            fprintf('Number of iterations completed: %d \n \n', i);
            break;
        end
    end
    
    if i == pm.outer_iter
        fprintf('Number of iterations completed: %d \n \n', i);
    end
end