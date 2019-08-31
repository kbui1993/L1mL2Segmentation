%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs PDHG to solve for u for each iteration of
%the DCA algorithm
%
%Input:
%   u0: initial iterate of the solution u
%   r: linear constant function
%   alpha: paramter for L1-alpha*L2
%   lambda: parameter
%   c: convexity parameter
%   inner_iter: number of iterations to run split Bregman
%   tau: step size for gradient descent of primal variable
%   sigma: step size for gradient ascent of dual variable
%
%Output:
%   u: solution of the DCA iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = PDHG(u0,r, alpha, lambda, c, inner_iter, tau, sigma)
    %obtain size
    [rows, cols] = size(u0);
    
    %set u0 to be the first iterate of solution
    u = u0;
    
    %obtain Dx and Dy of u0
    ux = Dx(u0);
    uy  =Dy(u0);
    
    %compute gradient of u
    ugrad = sqrt(abs(ux).^2+abs(uy).^2);
    eps = 1e-8;
    
    %compute qx and qy
    qx = ux./(ugrad+eps);
    qy = uy./(ugrad+eps);
    
    %preinitialize
    px = zeros(rows,cols);
    py = zeros(rows,cols);
    
    for i=1:inner_iter
        %store old u
        u_old = u;
        
        %update u
        u_new = (tau/(2*c*tau+1))*((2*c*u0+u/tau) - lambda*r+alpha*(Dxt(qx)+Dyt(qy))-Dxt(px)-Dyt(py));
        u_new = min(max(u_new,0),1);
        u_bar= 2*u_new - u;
        u = u_new;
        
        %update px and py
        px = px+sigma*Dx(u_bar);
        px = px./(max(abs(px),1));
        py = py+sigma*Dy(u_bar);
        py = py./(max(abs(py),1));
        
        % stop conditions
        relerr = norm(u_old - u)/max([norm(u_old), norm(u), eps]);
        
        if relerr < 1e-4 && i > 2
            fprintf('Number of iterations completed: %d \n \n', i);
            break;
        end
        
    end
end