fui8 = imread('shape.jpg');
f = double(fui8);
[N,M,~] = size(f);
fg = rescale_color_image(f);
fg(:,:,1) = add_noise2(fg(:,:,1), 0.1);
fg(:,:,2) = add_noise2(fg(:,:,2), 0.1);
fg(:,:,3) = add_noise2(fg(:,:,3), 0.1);


pm.outer_iter = 20;
pm.alpha = 0;
pm.lambda = 1;
pm.c = 0;
pm.inner_iter = 300;
pm.beta = 50;
pm.tau = 1/8;
pm.sigma = 1/8;
pm.method = 'PDHG';

u = make_circle(M,N,10);
u = double(u);

tic;
u1 = L1L2_color_two_phase(fg, u, pm);
%u1 = isoTV_color_two_phase(fg, u, pm);
time = toc
fg1 = fg;
fg1(:,:,1) = fg1(:,:,1).*double(u1>0.5);
fg1(:,:,2) = fg1(:,:,2).*double(u1>0.5);
fg1(:,:,3) = fg1(:,:,3).*double(u1>0.5);
figure; imagesc(fg1);