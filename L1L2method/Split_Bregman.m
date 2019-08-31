%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs split Bregman to solve for u for each iteration of
%the DCA algorithm
%
%Input:
%   u0: initial iterate of the solution u
%   r: linear constant function
%   alpha: paramter for L1-alpha*L2
%   lambda: parameter
%   c: convexity parameter
%   inner_iter: number of iterations to run split Bregman
%   beta: penalty parameter
%
%Output:
%   u: solution of the DCA iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = Split_Bregman(u0,r, alpha, lambda, c, inner_iter, beta)
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
    
    %compute kernel
    gg = fspecial('laplacian',0);
    uker = 2*c*eye(rows,cols)-beta*psf2otf(gg,[rows,cols]);
    
    %preinitialize
    tmp1 = zeros(size(u0));
    dx = zeros(rows,cols);
    dy = zeros(rows,cols);
    bx = zeros(rows,cols);
    by = zeros(rows,cols);
    for i = 1:inner_iter
        %store old u
        uold = u;
        
        %update u
        %rhs = 2*c*u0-lambda*r+beta*Dxt(dx-bx)+beta*Dyt(dy-by);
        %u = real(ifft2(ifftshift(fftshift(fft2(rhs))./uker)));
        tmp1(2:rows-1,2:cols-1) = u(3:rows,2:cols-1)+u(1:rows-2,2:cols-1)+u(2:rows-1,3:cols)+u(2:rows-1,1:cols-2);
        u = (beta/(2*c+4*beta))*(tmp1+Dx(dy-by)+Dy(dx-bx))+(2*c*u0)/(2*c+4*beta)-lambda*r/(2*c+4*beta);
        u = min(max(u,0),1);
        
        %compute Dx(u) and Dy(u)
        Dxu = Dx(u);
        Dyu = Dy(u);
        
        %update dx and dy
        dx = soft_threshold(Dxu+bx+alpha*qx/beta, 1/beta);
        dy = soft_threshold(Dyu+by+alpha*qy/beta, 1/beta);
        
        %update bx and by
        bx = bx+Dxu-dx;
        by = by+Dyu-dy;
        
        % stop conditions
        relerr = norm(uold - u)/max([norm(uold), norm(u), eps]);
        
        if relerr < 1e-4 && i > 2
            fprintf('Number of iterations completed: %d \n \n', i);
            break;
        end
    end
    
end