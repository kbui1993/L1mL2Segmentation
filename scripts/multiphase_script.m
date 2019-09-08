fui8 = imread('mri2.jpg');
f = double(fui8);
[N,M] = size(f);
fnoise = add_noise2(f,0);
fg = double(fnoise);
fg = rescale_image(fg);


pm.c1 = 0.4;
pm.c2 = 0.8;
pm.c3 = 0.9;
pm.c4 = 0.;
pm.outer_iter = 5;
pm.alpha = 0.1;
pm.lambda = 100;
pm.c = 1e-4;
pm.inner_iter = 300;
pm.beta = 50;
pm.tau = 1/6;
pm.sigma = 1/6;
pm.method = 'PDHG';

u1 = make_circle_shift_x(M,N,10, -5);
u2 = make_circle_shift_x(M,N,10, 5);
u1 = double(u1);
u2 = double(u2);
% u1 = zeros(M,N);
% u2 = zeros(M,N);
% for i=1:M
%     for j=1:N
%         u1(i,j) = sin((pi*i)/3)*sin((pi*j)/3);
%         u2(i,j) = sin((pi*i)/10)*sin((pi*j)/10);
%     end
% end

tic;
[U1,U2] = L1L2_four_phase(fg, u1, u2, pm);
time = toc
figure;
subplot(2,3,1); imagesc(f); axis off; axis square; colormap gray; title('Original');
subplot(2,3,2); imagesc(double(U1>0.5).*double(U2>0.5)); axis off; axis square; title('Phase 1');
subplot(2,3,3); imagesc(double(U1>0.5).*double(U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(2,3,5); imagesc(double(U1<=0.5).*double(U2>0.5)); axis off; axis square; title('Phase 3');
subplot(2,3,6); imagesc(double(U1<=0.5).*double(U2<=0.5)); axis off; axis square; title('Phase 4');