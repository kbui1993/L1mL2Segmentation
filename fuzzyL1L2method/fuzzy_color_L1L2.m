%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs L1-alpha* L2 TV for fuzzy region competition model
%for color segmentation.
%Input:
%   f: image
%   u_initial: initialization
%   pm: set of parameters
%   k: number of regions
%
%Output:
%   u_result: segmentation result
%   c_result: color values for each region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_result, c_result] = fuzzy_color_L1L2(f, u_initial, pm, k)

%separate the channels of f
f_r = f(:,:,1);
f_g = f(:,:,2);
f_b = f(:,:,3);
    

%obtain dimensions of the image
M = size(f,1);
N = size(f,2);

%preinitialize cell
u_result = u_initial;

%sum up the u's
sum_u = sum(cat(3,u_result{:}),3);

%normalize each u
for i = 1:k
    u_result{i} = u_result{i}./sum_u;
end

%set error
error = 10^(-10)*eye(M,N);

%compute constant vectors for each region
c_result = {};
rel_err = ones(k,1);
for i = 1:k
    c_result{i}=[sum(f_r(:).*u_result{i}(:))./sum(u_result{i}(:)+error(:));
        sum(f_g(:).*u_result{i}(:))./sum(u_result{i}(:)+error(:));
        sum(f_b(:).*u_result{i}(:))./sum(u_result{i}(:)+error(:))];
end

%optimize model
for i =1:pm.outer_iter
    
    %set old_u_result to be current u_result
    old_u_result = u_result;
    
    %minimize for each u_j
    for j = 1:k
        
        %get the constant vector for region j
        c = c_result{j};
        
        %compute the fidelity function value
        f_k = (f_r-c(1)).^2 + (f_g-c(2)).^2+(f_b-c(3)).^2;
        
        %compute 1-sum_{i \neq j} u_i
        u2 = ones(M,N) - (sum_u-u_result{j});
        
        %perform primal-dual hybrid gradient with line search
        u_result{j} = PDHGLS2(u_result{j}, u2, f_k, pm.alpha, pm.lambda, pm.nu, pm.c, pm.inner_iter, pm.tau, pm.beta);
        
        %recompute sum u_i
        sum_u = sum(cat(3,u_result{:}),3);
    end
    
    %update constant vectors for each region
    for j = 1:k
        c_result{j}=[sum(f_r(:).*u_result{j}(:))./sum(u_result{j}(:)+error(:));
            sum(f_g(:).*u_result{j}(:))./sum(u_result{j}(:)+error(:));
            sum(f_b(:).*u_result{j}(:))./sum(u_result{j}(:)+error(:))];
        rel_err(j) = norm(old_u_result{j} - u_result{j}, 'fro')/max([norm(old_u_result{j},'fro'), norm(u_result{j}, 'fro'), eps]);
    end
    
    %check stopping condition
    if max(rel_err) < 1e-4 && i >2
        fprintf('Number of DCA inner iterations completed: %d \n \n', i);
        break;
    end
end

%print number of outer iterations
if i == pm.outer_iter
    fprintf('Number of DCA outer iterations completed: %d \n \n', i);
end

end