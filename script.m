fui8 = imread('cameraman.tif');
f = double(fui8);
[N,M] = size(f);
fnoise = add_noise2(f,0.1);
fg = double(fnoise);
fg = rescale_image(fg);


pm.c1 = 0.1;
pm.c2 = 0.75;
pm.outer_iter = 20;
pm.alpha = 0.1;
pm.lambda = 5;
pm.c = 1e-4;
pm.inner_iter = 300;
pm.beta = 100;
pm.tau = 1/6;
pm.sigma = 1/6;
pm.method = 'PDHG';

%u = make_circle(M,N,10);
%u = double(u);
u = zeros(N,M);
for i = 1:N
    for j=1:M
        u(i,j) = sin(pi*i/N)*sin(pi*j/M);
    end
end

tic;
u1 = L1L2_two_phase(fg, u, pm);
time = toc
figure; imagesc(double(u1>0));