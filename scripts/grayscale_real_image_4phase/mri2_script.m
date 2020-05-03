%loading mri2.jpg image
fui8 = imread('mri2.jpg');
f = double(fui8);

%rescale image and get image size
fg = rescale_image(f);
[N,M] = size(fg);

%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 10000;
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
[L1L2_U1,L1L2_U2,L1L2_c1,L1L2_c2,L1L2_c3,L1L2_c4] = L1L2_four_phase(fg, u1, u2, pm);
toc
L1L2_approx_im = double(L1L2_U1>0.5).*double(L1L2_U2>0.5).*L1L2_c1 + double(L1L2_U1>0.5).*double(L1L2_U2<=0.5).*L1L2_c2 +...
    double(L1L2_U1<=0.5).*double(L1L2_U2>0.5).*L1L2_c3 + double(L1L2_U1<=0.5).*double(L1L2_U2<=0.5).*L1L2_c4;

%L1-0.5L2 segmentation
pm.alpha = 0.5;
tic;
[L1L2_05_U1,L1L2_05_U2,L1L2_05_c1,L1L2_05_c2,L1L2_05_c3,L1L2_05_c4] = L1L2_four_phase(fg, u1, u2, pm);
toc
L1L2_05_approx_im = double(L1L2_05_U1>0.5).*double(L1L2_05_U2>0.5).*L1L2_05_c1 + double(L1L2_05_U1>0.5).*double(L1L2_05_U2<=0.5).*L1L2_05_c2 +...
    double(L1L2_05_U1<=0.5).*double(L1L2_05_U2>0.5).*L1L2_05_c3 + double(L1L2_05_U1<=0.5).*double(L1L2_05_U2<=0.5).*L1L2_05_c4;

%anisotropic CV
pm.alpha = 0;
pm.c = 0;
tic;
[ani_U1,ani_U2, ani_c1, ani_c2, ani_c3, ani_c4] = L1L2_four_phase(fg, u1, u2, pm);
toc

ani_approx_im = double(ani_U1>0.5).*double(ani_U2>0.5).*ani_c1 + double(ani_U1>0.5).*double(ani_U2<=0.5).*ani_c2 +...
    double(ani_U1<=0.5).*double(ani_U2>0.5).*ani_c3 + double(ani_U1<=0.5).*double(ani_U2<=0.5).*ani_c4;


%isotropic CV
tic;
[iso_U1,iso_U2, iso_c1, iso_c2, iso_c3, iso_c4] = isoTV_four_phase(fg, u1, u2, pm);
toc

iso_approx_im = double(iso_U1>0.5).*double(iso_U2>0.5).*iso_c1 + double(iso_U1>0.5).*double(iso_U2<=0.5).*iso_c2 +...
    double(iso_U1<=0.5).*double(iso_U2>0.5).*iso_c3 + double(iso_U1<=0.5).*double(iso_U2<=0.5).*iso_c4;


%plot figure
figure;
subplot(4,6,1); imagesc(fg); axis off; axis square; colormap gray; title('Original');
subplot(4,6,2); imagesc(double(L1L2_U1>0.5).*double(L1L2_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,6,3); imagesc(double(L1L2_U1>0.5).*double(L1L2_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,6,4); imagesc(double(L1L2_U1<=0.5).*double(L1L2_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,6,5); imagesc(double(L1L2_U1<=0.5).*double(L1L2_U2<=0.5)); axis off; axis square; title('Phase 4');
subplot(4,6,6); imagesc(L1L2_approx_im); axis off; axis square; title('Approximation')

subplot(4,6,8); imagesc(double(L1L2_05_U1>0.5).*double(L1L2_05_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,6,9); imagesc(double(L1L2_05_U1>0.5).*double(L1L2_05_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,6,10); imagesc(double(L1L2_05_U1<=0.5).*double(L1L2_05_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,6,11); imagesc(double(L1L2_05_U1<=0.5).*double(L1L2_05_U2<=0.5)); axis off; axis square; title('Phase 4');
subplot(4,6,12); imagesc(L1L2_05_approx_im); axis off; axis square; title('Approximation')

subplot(4,6,14); imagesc(double(ani_U1>0.5).*double(ani_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,6,15); imagesc(double(ani_U1>0.5).*double(ani_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,6,16); imagesc(double(ani_U1<=0.5).*double(ani_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,6,17); imagesc(double(ani_U1<=0.5).*double(ani_U2<=0.5)); axis off; axis square; title('Phase 4');
subplot(4,6,18); imagesc(ani_approx_im); axis off; axis square; title('Approximation')

subplot(4,6,20); imagesc(double(iso_U1>0.5).*double(iso_U2>0.5)); axis off; axis square; title('Phase 1');
subplot(4,6,21); imagesc(double(iso_U1>0.5).*double(iso_U2<=0.5)); axis off; axis square; title('Phase 2');
subplot(4,6,22); imagesc(double(iso_U1<=0.5).*double(iso_U2>0.5)); axis off; axis square; title('Phase 3');
subplot(4,6,23); imagesc(double(iso_U1<=0.5).*double(iso_U2<=0.5)); axis off; axis square; title('Phase 4');
subplot(4,6,24); imagesc(iso_approx_im); axis off; axis square; title('Approximation')
