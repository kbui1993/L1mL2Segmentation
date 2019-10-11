%fui8 = imread('cameraman.tif');
fui8 = imread('angio.jpg');
fui8 = rgb2gray(fui8);
f = double(fui8);
[N,M] = size(f);
fnoise = add_noise2(f,0);
fg = double(fnoise);
fg = rescale_image(fg);

pm.outer_iter = 2;
pm.alpha = 0;
pm.lambda = 10;
pm.c = 0;
pm.inner_iter = 300;
pm.beta = 100;
pm.tau = 1/8;
pm.sigma = 1/8;
pm.method = 'PDHG';

u = make_circle(M,N,10);
u = double(u);
% u = zeros(N,M);
% for i = 1:N
%     for j=1:M
%         u(i,j) = sin(pi*i/N)*sin(pi*j/M);
%     end
% end

tic;
u1 = L1L2_two_phase(fg, u, pm);
%u1 = isoTV_two_phase(fg, u, pm);
time = toc
figure; imagesc(double(u1>0.5));