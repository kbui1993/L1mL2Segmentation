%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs primal-dual hybrid gradient with line search for
%the fuzzy region competition model.
%Input:
%   u0: initialization for the function u_j that needs to be minimized
%   u0_2: 1- sum_{i \neq j} u_i
%   r: linear constant function associated with u_j
%   alpha: parameter for L1-\alpha*L2
%   lambda: fidelity parameter
%   nu: quadratic penalty parameter for sum-of-one constraint
%   c: strong convexity parameter
%   inner_iter: number of max iterations to perform PDHGLS
%   tau: initial step-size for the primal descent
%   beta: primal-dual step size ratio
%
%Output:
%   u: solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = PDHGLS2(u0, u0_2, r, alpha, lambda, nu, c, inner_iter, tau, beta)
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
        u_new = (tau_old/(2*c*tau_old+1+nu*tau_old))*((2*c*u0+u/tau_old+nu*u0_2)- lambda*r+alpha*(Dxt(qx)+Dyt(qy))-Dxt(px)-Dyt(py));
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
        
        % stopping conditions
        relerr = norm(u_old - u, 'fro')/max([norm(u_old, 'fro'), norm(u, 'fro'), eps]);
        
        if relerr < 1e-6 && i > 2
            fprintf('Number of PDHGLS iterations completed: %d \n \n', i);
            break;
        end
        
    end
end