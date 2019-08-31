function u = L1L2_two_phase(f,u_initial, pm)
    %compute r
    r = (pm.c1-f).^2- (pm.c2-f).^2;
    
    %set u
    u = u_initial;
    
    %run DCA algorithm
    for i = 1:pm.outer_iter
        u_old = u;
        if strcmp(pm.method, 'PDHG')
            u = PDHG(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.sigma);
        elseif strcmp(pm.method, 'aPDHG')
            u = aPDHG(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.tau, pm.sigma);
        else
            u = Split_Bregman(u, r, pm.alpha, pm.lambda, pm.c, pm.inner_iter, pm.beta);
        end
        relerr = norm(u_old - u)/max([norm(u_old), norm(u), eps]);
        if relerr <5e-3 && i > 2
            fprintf('Number of iterations completed: %d \n \n', i);
            break;
        end
    end
    
    if i == pm.outer_iter
        fprintf('Number of iterations completed: %d \n \n', i);
    end
end