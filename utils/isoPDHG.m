%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs PDHG to solve for u for each iteration of
%the DCA algorithm
%
%Input:
%   u0: initial iterate of the solution u
%   r: linear constant function
%   lambda: parameter
%   inner_iter: number of iterations to run PDHG
%   tau: step size for gradient descent of primal variable
%   sigma: step size for gradient ascent of dual variable
%
%Output:
%   u: solution of the DCA iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = isoPDHG(u0,r, lambda, inner_iter, tau, sigma)
    %obtain size
    [rows, cols] = size(u0);
    
    %set u0 to be the first iterate of solution
    u = u0;
    eps = 1e-8;
    
    %preinitialize
    px = zeros(rows,cols);
    py = zeros(rows,cols);
    
    for i=1:inner_iter
        %store old u
        u_old = u;
        
        %update u
        u_new = u - tau*(lambda*r+Dxt(px)+Dyt(py));
        u_new = min(max(u_new,0),1);
        u_bar= 2*u_new - u;
        u = u_new;
        
        %update px and py
        px = px+sigma*Dx(u_bar);
        py = py+sigma*Dy(u_bar);
        norm1 = sqrt(max(px.*px+py.*py,1));
        px = px./norm1;
        py = py./norm1;
        
        % stop conditions
        relerr = norm(u_old - u)/max([norm(u_old), norm(u), eps]);
        
        if relerr < 1e-6 && i > 2
            fprintf('Number of iterations completed: %d \n \n', i);
            break;
        end
        
    end
end