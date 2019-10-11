%loading in angio image and rescaling it to [0,1]
fui8 = imread('shape.png');
fui8 = rgb2gray(fui8);
f = double(fui8);
[N,M] = size(f);
% fnoise = add_noise2(f,0.5);
% fg = double(fnoise);
% fg = rescale_image(fg);

fblur = imgaussfilt(f,2);
fg = rescale_image(fblur);

f_clean = rescale_image(f);

%setting up parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 50;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.sigma = 1/8;
pm.method = 'PDHG';

%preinitialize u
u = make_circle(M,N,10);
u = double(u);

tic;
L1L2_1_u1 = L1L2_two_phase(fg, u, pm);
time = toc

pm.alpha = 0.5;
tic;
L1L2_pt5_u1 = L1L2_two_phase(fg, u, pm);
time = toc

pm.alpha = 0;
pm.c = 0;
tic;
ani_u1 = L1L2_two_phase(fg, u, pm);
time = toc

tic;
iso_u1 = isoTV_two_phase(fg, u, pm);
time = toc

figure;
subplot(2,3,1); imagesc(fg); axis off; axis square; colormap gray; title('Original');
subplot(2,3,2); imagesc(double(L1L2_1_u1>0.5)); axis off; axis square; title('L1-L2');
subplot(2,3,3); imagesc(double(L1L2_pt5_u1>0.5)); axis off; axis square; title('L1-0.5L2');
subplot(2,3,5); imagesc(double(ani_u1>0.5)); axis off; axis square; title('anisotropic');
subplot(2,3,6); imagesc(double(iso_u1>0.5)); axis off; axis square; title('isotropic')

a1 = double(L1L2_1_u1>0.5);
a2 = double(L1L2_pt5_u1>0.5);
a3 = double(ani_u1>0.5);
a4 = double(iso_u1>0.5);

ssim(a1, f_clean)
ssim(a2, f_clean)
ssim(a3, f_clean)
ssim(a4, f_clean)