%This script runs experiments to compare the proposed two-phase method 
%against isotropic and anisotropic Chan Vese segmentation on synthetic
%RGB image with Gaussian blur.

%generate color synthetic image
color_synthetic_image;
f = M;
foriginal=M;

%get image size
[N,M,~] = size(f);

%apply blur to each channel
f(:,:,1) = imgaussfilt(f(:,:,1), 2);
f(:,:,2) = imgaussfilt(f(:,:,2), 2);
f(:,:,3) = add_noise2(f(:,:,3), 2);
fg = rescale_color_image(f);

%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 300;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.sigma = 1/8;
pm.method = 'PDHG';

%set image segmentation initialization
u = make_circle(M,N,10);
u = double(u);

%L1-L2
tic;
L1_L2_u1 = L1L2_color_two_phase(fg, u, pm);
time = toc

%L1-0.5L2
pm.alpha =0.5;
tic;
L1_0pt5_L2_u1 = L1L2_color_two_phase(fg, u, pm);
time = toc

%anisotropic CV
pm.alpha = 0;
pm.c = 0;
tic;
L1_u1 = L1L2_color_two_phase(fg, u, pm);
time = toc

%isotropic CV
tic;
iso_u1 = isoTV_color_two_phase(fg, u, pm);
time = toc

%%reconstruct image from image segmentation

%L1-L2
L1_L2M1=ones(n,n);
L1_L2M1(L1_L2_u1<0.5)=0;
L1_L2M2=double(L1_L2_u1<0.5)*90;
L1_L2M3=double(L1_L2_u1<0.5)*150;
L1_L2M = zeros(n,n,3);
L1_L2M(:,:,1)=L1_L2M1;
L1_L2M(:,:,2) = L1_L2M2;
L1_L2M(:,:,3) = L1_L2M3;

%L1-0.5L2
L1_0pt5_L2M1=ones(n,n);
L1_0pt5_L2M1(L1_0pt5_L2_u1<0.5)=0;
L1_0pt5_L2M2=double(L1_0pt5_L2_u1<0.5)*90;
L1_0pt5_L2M3=double(L1_0pt5_L2_u1<0.5)*150;
L1_0pt5_L2M = zeros(n,n,3);
L1_0pt5_L2M(:,:,1)=L1_0pt5_L2M1;
L1_0pt5_L2M(:,:,2) = L1_0pt5_L2M2;
L1_0pt5_L2M(:,:,3) = L1_0pt5_L2M3;

%anisotropic
L1M1=ones(n,n);
L1M1(L1_u1<0.5)=0;
L1M2=double(L1_u1<0.5)*90;
L1M3=double(L1_u1<0.5)*150;
L1M = zeros(n,n,3);
L1M(:,:,1)=L1M1;
L1M(:,:,2) = L1M2;
L1M(:,:,3) = L1M3;

%isotropic
isoM1=ones(n,n);
isoM1(iso_u1<0.5)=0;
isoM2=double(iso_u1<0.5)*90;
isoM3=double(iso_u1<0.5)*150;
isoM = zeros(n,n,3);
isoM(:,:,1)=isoM1;
isoM(:,:,2) = isoM2;
isoM(:,:,3) = isoM3;

%compute ssim
ssim(L1_L2M,foriginal)
ssim(L1_0pt5_L2M, foriginal)
ssim(L1M, foriginal)
ssim(isoM, foriginal)

%plot segmentation
figure;
subplot(2,3,1); imagesc(f); axis off; axis square; title('Original');
subplot(2,3,2); imagesc(f); hold on; contour(double(L1_L2_u1>0.5), 'k'); axis off; axis square; title('L1-L2');
subplot(2,3,3); imagesc(f); hold on; contour(double(L1_0pt5_L2_u1>0.5), 'k'); axis off; axis square; title('L1-0.5L2');
subplot(2,3,5); imagesc(f); hold on; contour(double(L1_u1>0.5), 'k'); axis off; axis square; title('Anisotropic');
subplot(2,3,6); imagesc(f); hold on; contour(double(iso_u1>0.5), 'k'); axis off; axis square; title('Isotropic');
