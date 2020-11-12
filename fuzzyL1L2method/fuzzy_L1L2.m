%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs L1-alpha* L2 TV for fuzzy region competition model
%for grayscale segmentation.
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
function [u_result, c_result] = fuzzy_L1L2(f, u_initial, pm, k)

%obtain dimensions of the image
M = size(f,1);
N = size(f,2);

%preinitialize cell
u_result = u_initial;

%sum up the u's
sum_u = sum(cat(3,u_result{:}),3);

%normalize each u's
for i = 1:k
    u_result{i} = u_result{i}./sum_u;
end

%set error
error = 10^(-10)*eye(M,N);

%compute constant values for each region
c_result = {};
rel_err = ones(k,1);
for i = 1:k
    c_result{i} = sum(f(:).*u_result{i}(:))./sum(u_result{i}(:) + error(:));
end

%optimize model
for i =1:pm.outer_iter
    
    %set old u_result to be current u result
    old_u_result = u_result;
    
    %perform inner convex optimization for each u_j
    for j = 1:k
        
        %compute fidelity for region j
        f_k = (f-c_result{j}).^2;
        
        %compute 1-sum_{i \neq j} u_i
        u2 = ones(M,N) - (sum_u-u_result{j});
        
        %perform primal-dual hybrid gradient with line search
        u_result{j} = PDHGLS2(u_result{j}, u2, f_k, pm.alpha, pm.lambda, pm.nu, pm.c, pm.inner_iter, pm.tau, pm.beta);
        
        %compute new sum across u_j's
        sum_u = sum(cat(3,u_result{:}),3);
    end
    
    %update the c_j's
    for j = 1:k
        c_result{j} = sum(f(:).*u_result{j}(:))./sum(u_result{j}(:) + error(:));
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