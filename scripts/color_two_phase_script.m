fui8 = imread('113044.jpg');
f = double(fui8);
[N,M,~] = size(f);
fg = rescale_color_image(f);

pm.outer_iter = 20;
pm.alpha = 0.1;
pm.lambda = 1000;
pm.c = 1e-4;
pm.inner_iter = 300;
pm.beta = 50;
pm.tau = 1/4;
pm.sigma = 1/4;
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
u1 = L1L2_color_two_phase(fg, u, pm);
time = toc
fg1 = fg;
fg1(:,:,1) = fg1(:,:,1).*double(u1>0.5);
fg1(:,:,2) = fg1(:,:,2).*double(u1>0.5);
fg1(:,:,3) = fg1(:,:,3).*double(u1>0.5);
figure; imagesc(fg1);