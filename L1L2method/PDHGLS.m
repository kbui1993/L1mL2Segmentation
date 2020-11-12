%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs PDHG with line search to solve for u for each iteration of
%the DCA algorithm
%
%Input:
%   u0: initial iterate of the solution u
%   r: linear constant function
%   alpha: paramter for L1-alpha*L2
%   lambda: fidelity parameter
%   c: strong convexity parameter
%   inner_iter: number of iterations to run PDHG
%   tau: step size for gradient descent of primal variable
%   beta: ratio parameter
%
%Output:
%   u: solution of the DCA iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = PDHGLS(u0,r, alpha, lambda, c, inner_iter, tau, beta)
    %set decreasing parameter for line search
    mu = 7.5e-5;
    
    %set delta parameter for when to end line search
    delta = 0.9999;
    
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
    
    %initialize changing step size
    tau_old = tau;
    theta = 1;
    for i=1:inner_iter
        %store old u
        u_old = u;
        
        %update u
        u_new = (tau_old/(2*c*tau_old+1))*((2*c*u0+u/tau_old) - lambda*r+alpha*(Dxt(qx)+Dyt(qy))-Dxt(px)-Dyt(py));
        u_new = min(max(u_new,0),1);
        u = u_new;
        
        %performing line search starts here
        tau_new = sqrt(1+theta)*tau_old;
        theta = tau_new/tau_old;
        
        u_bar= (1+theta)*u_new-theta*u_old;
        
        %update px and py
        px_new = px+beta*tau_new*Dx(u_bar);
        px_new = px_new./(max(abs(px_new),1));
        py_new = py+beta*tau_new*Dy(u_bar);
        py_new = py_new./(max(abs(py_new),1));
        
        while(sqrt(beta)*tau_new*(norm([Dxt(px_new);Dyt(py_new)]-[Dxt(px);Dyt(py)],'fro'))>delta*(norm([px_new;py_new]-[px;py],'fro')))
            %update tau
            tau_new = tau_new*mu;
            
            %compute thetha
            theta = tau_new/tau_old;
            
            %compute u_bar
            u_bar= (1+theta)*u_new-theta*u_old;
            
            %update px and py
            px_new = px+beta*tau_new*Dx(u_bar);
            px_new = px_new./(max(abs(px_new),1));
            py_new = py+beta*tau_new*Dy(u_bar);
            py_new = py_new./(max(abs(py_new),1));
        end
        
        %update
        tau_old = tau_new;
        px = px_new;
        py = py_new;
        
        % stop conditions
        relerr = norm(u_old - u, 'fro')/max([norm(u_old, 'fro'), norm(u, 'fro'), eps]);
        
        if relerr < 1e-6 && i > 2
            fprintf('Number of PDHGLS iterations completed: %d \n \n', i);
            break;
        end
        
    end
end