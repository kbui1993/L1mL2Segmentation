%loading hand.jpg image
fui8 = imread('hand.jpg');

%convert to grayscale image
fui8 = rgb2gray(fui8);
f = double(fui8);

%rescale image and get image size
fg = rescale_image(f);
[N,M] = size(fg);

%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 50;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/4;
pm.sigma = 1/4;
pm.method = 'PDHG';

%set segmentation initialization
u1 = make_circle_shift_x(M,N,10, -5);
u2 = make_circle_shift_x(M,N,10, 5);
u1 = double(u1);
u2 = double(u2);

%L1-L2 segmentation
tic;
[L1L2_U1,L1L2_U2] = L1L2_four_phase(fg, u1, u2, pm);
toc

%L1-0.5L2 segmentation
pm.alpha = 0.5;
tic;
[L1L2_05_U1,L1L2_05_U2] = L1L2_four_phase(fg, u1, u2, pm);
toc

%anisotropic CV
pm.alpha = 0;
pm.c = 0;
tic;
[ani_U1,ani_U2] = L1L2_four_phase(fg, u1, u2, pm);
toc

%isotropic CV
tic;
[iso_U1,iso_U2] = isoTV_four_phase(fg, u1, u2, pm);
toc

%plot figure
figure;
subplot(4,5,1); imagesc(fg); axis off; axis square; colormap gray; title('Original');
subplot(4,5,2); imagesc(double(L1L2_U1>0.5).*double(L1L2_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,5,3); imagesc(double(L1L2_U1>0.5).*double(L1L2_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,5,4); imagesc(double(L1L2_U1<=0.5).*double(L1L2_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,5,5); imagesc(double(L1L2_U1<=0.5).*double(L1L2_U2<=0.5)); axis off; axis square; title('Phase 4');

subplot(4,5,7); imagesc(double(L1L2_05_U1>0.5).*double(L1L2_05_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,5,8); imagesc(double(L1L2_05_U1>0.5).*double(L1L2_05_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,5,9); imagesc(double(L1L2_05_U1<=0.5).*double(L1L2_05_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,5,10); imagesc(double(L1L2_05_U1<=0.5).*double(L1L2_05_U2<=0.5)); axis off; axis square; title('Phase 4');

subplot(4,5,12); imagesc(double(ani_U1>0.5).*double(ani_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,5,13); imagesc(double(ani_U1>0.5).*double(ani_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,5,14); imagesc(double(ani_U1<=0.5).*double(ani_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,5,15); imagesc(double(ani_U1<=0.5).*double(ani_U2<=0.5)); axis off; axis square; title('Phase 4');

subplot(4,5,17); imagesc(double(iso_U1>0.5).*double(iso_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,5,18); imagesc(double(iso_U1>0.5).*double(iso_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,5,19); imagesc(double(iso_U1<=0.5).*double(iso_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,5,20); imagesc(double(iso_U1<=0.5).*double(iso_U2<=0.5)); axis off; axis square; title('Phase 4');