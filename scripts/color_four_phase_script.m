fui8 = imread('color_shape.jpg');
f = double(fui8);
[N,M,~] = size(f);
f(:,:,1) = add_noise2(f(:,:,1), 0.2);
f(:,:,2) = add_noise2(f(:,:,2), 0.2);
f(:,:,3) = add_noise2(f(:,:,3), 0.2);
fg = rescale_color_image(f);
fg = double(fg);

pm.outer_iter = 20;
pm.alpha = 0.1;
pm.lambda = 10;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.sigma = 1/8;
pm.method = 'PDHG';

u1 = make_circle_shift_x(M,N,10, -5);
u2 = make_circle_shift_x(M,N,10, 5);
u1 = double(u1);
u2 = double(u2);
% u = zeros(N,M);
% for i = 1:N
%     for j=1:M
%         u(i,j) = sin(pi*i/N)*sin(pi*j/M);
%     end
% end

tic;
[U1, U2] = L1L2_color_four_phase(fg, u1, u2, pm);
%[U1, U2] = isoTV_color_four_phase(fg, u1, u2, pm);
time = toc

a1 = double(U1>0.5).*double(U2>0.5);
fg1 = fg;
fg1(:,:,1) = fg1(:,:,1).*a1;
fg1(:,:,2) = fg1(:,:,2).*a1;
fg1(:,:,3) = fg1(:,:,3).*a1;

a2 = double(U1>0.5).*double(U2<=0.5);
fg2 = fg;
fg2(:,:,1) = fg2(:,:,1).*a2;
fg2(:,:,2) = fg2(:,:,2).*a2;
fg2(:,:,3) = fg2(:,:,3).*a2;

a3 = double(U1<=0.5).*double(U2>0.5);
fg3 = fg;
fg3(:,:,1) = fg3(:,:,1).*a3;
fg3(:,:,2) = fg3(:,:,2).*a3;
fg3(:,:,3) = fg3(:,:,3).*a3;

a4 = double(U1<=0.5).*double(U2<=0.5);
fg4 = fg;
fg4(:,:,1) = fg4(:,:,1).*a4;
fg4(:,:,2) = fg4(:,:,2).*a4;
fg4(:,:,3) = fg4(:,:,3).*a4;
figure;
subplot(2,3,1); imagesc(fg); axis off; axis square; title('Original');
subplot(2,3,2); imagesc(fg1); axis off; axis square; title('Phase 1');
subplot(2,3,3); imagesc(fg2); axis off; axis square; title('Phase 2');
subplot(2,3,5); imagesc(fg3); axis off; axis square; title('Phase 3');
subplot(2,3,6); imagesc(fg4); axis off; axis square; title('Phase 4');