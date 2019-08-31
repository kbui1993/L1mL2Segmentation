fui8 = imread('cameraman.tif');
f = double(fui8);
[N,M] = size(f);
fnoise = add_noise2(f,0);
fg = double(fnoise);
fg = rescale_image(fg);


pm.c1 = 0.1;
pm.c2 = 0.75;
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 5;
pm.c = 1e-9;
pm.inner_iter = 100;
pm.beta = 5;
pm.tau = 1/4;
pm.sigma = 1/4;
pm.method = 'PDHG';

u = make_circle(M,N,50);
u = double(u);

tic;
u1 = two_phase(fg, u, pm);
time = toc