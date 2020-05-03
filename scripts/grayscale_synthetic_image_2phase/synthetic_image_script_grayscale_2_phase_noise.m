%This script runs experiments to compare the proposed two-phase method 
%against isotropic and anisotropic Chan Vese segmentation on synthetic
%grayscale image with Gaussian noise.

%read image
fui8 = imread('shape.png');
fui8 = rgb2gray(fui8);
f = double(fui8);
[N,M] = size(f);

%%with Gaussian noise
rng(1234);
fnoise = add_noise2(f,0.20);

%with salt and pepper noise
%fnoise = imnoise(f,'salt & pepper', 0.20);
fg = double(fnoise);
fg = rescale_image(fg);

%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 2.5;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.sigma = 1/8;
pm.method = 'PDHG';

%preinitialize segmentation
u = make_circle(M,N,10);

u = double(u);


%L1-1.0L2
tic;
L1_L2_u1 = L1L2_two_phase(fg, u, pm);
time = toc

%L1 - 0.5L2
pm.alpha = 0.5;
tic;
L1_0pt5_L2_u1 = L1L2_two_phase(fg, u, pm);
toc

%anisotropic CV
pm.alpha =0;
pm.c = 0;
tic;
L1_u1 = L1L2_two_phase(fg, u, pm);
toc

%isotropic CV
tic;
iso_u1 = isoTV_two_phase(fg, u, pm);
toc

fg1 = rescale_image(f);

%compute dice
dice(double(L1_L2_u1>0.5), fg1)
dice(double(L1_0pt5_L2_u1>0.5), fg1)
dice(double(L1_u1>0.5), fg1)
dice(double(iso_u1>0.5), fg1)

%plot segmentation
figure;
subplot(2,3,1); imagesc(fg); axis off; axis square; colormap gray; title('Original');
subplot(2,3,2); imagesc(fg); hold on; contour(double(L1_L2_u1>0.5), 'g'); axis off; axis square; title('L1-L2');
subplot(2,3,3); imagesc(fg); hold on; contour(double(L1_0pt5_L2_u1>0.5), 'g'); axis off; axis square; title('L1-0.5L2');
subplot(2,3,5); imagesc(fg); hold on; contour(double(L1_u1>0.5), 'g'); axis off; axis square; title('Anisotropic');
subplot(2,3,6); imagesc(fg); hold on; contour(double(iso_u1>0.5), 'g'); axis off; axis square; title('Isotropic');

