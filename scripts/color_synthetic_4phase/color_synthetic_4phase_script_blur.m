%This script runs experiments to compare the proposed four-phase method 
%against isotropic and anisotropic Chan Vese segmentation on synthetic
%RGB image with Gaussian blur.

%generate synthetic color image
color_synthetic_image2;
f = double(M);
[N,M,~] = size(f);

%save rescaled image
fg1 = rescale_color_image(f);

%apply Gaussian blur
f(:,:,1) = imgaussfilt(f(:,:,1), 3);
f(:,:,2) = imgaussfilt(f(:,:,2), 3);
f(:,:,3) = imgaussfilt(f(:,:,3), 3);
fg = rescale_color_image(f);
fg = double(fg);

%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda =5.425;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/6;
pm.sigma = 1/6;
pm.method = 'PDHG';

%set image segmentation initialization
u1 = make_circle_shift_x(M,N,10, -5);
u2 = make_circle_shift_x(M,N, 10, 5);
u1 = double(u1);
u2 = double(u2);

%L1-L2
tic;
[L1L2_U1,L1L2_U2] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

%L1-0.5L2
pm.alpha = 0.5;
tic;
[L1L2_05_U1,L1L2_05_U2] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

%anisotropic
pm.alpha = 0;
pm.c = 0;
tic;
[ani_U1,ani_U2] = L1L2_color_four_phase(fg, u1, u2, pm);
toc

%isotropic
tic;
[iso_U1,iso_U2] = isoTV_color_four_phase(fg, u1, u2, pm);
toc

%%reconstruct image from image segmentation
%L1-L2
L1_L2M1=ones(n,n);
L1_L2M3=zeros(n,n);
L1_L2M1(double(L1L2_U1>0.5).*double(L1L2_U2>0.5)==1 | double(L1L2_U1>0.5).*double(L1L2_U2<=0.5)==1)=0.7;
L1_L2M2 = (double(L1L2_U1>0.5).*double(L1L2_U2>0.5)==1)*0.1+ (double(L1L2_U1<=0.5).*double(L1L2_U2<=0.5)==1)*0.7;
L1_L2M3(double(L1L2_U1>0.5).*double(L1L2_U2<=0.5)==1 | double(L1L2_U1<=0.5).*double(L1L2_U2<=0.5)==1)=0.1;

L1_L2M = zeros(n,n,3);
L1_L2M(:,:,1)=L1_L2M1;
L1_L2M(:,:,2) = L1_L2M2;
L1_L2M(:,:,3) = L1_L2M3;

L1_L2M = rescale_color_image(L1_L2M);

%L1-0.5L2
L1_0pt5_L2M1=ones(n,n);
L1_0pt5_L2M3=zeros(n,n);
L1_0pt5_L2M1(double(L1L2_05_U1>0.5).*double(L1L2_05_U2>0.5)==1 | double(L1L2_05_U1>0.5).*double(L1L2_05_U2<=0.5)==1)=0.7;
L1_0pt5_L2M2 = (double(L1L2_05_U1>0.5).*double(L1L2_05_U2>0.5)==1)*0.1+ (double(L1L2_05_U1<=0.5).*double(L1L2_05_U2>0.5)==1)*0.7;
L1_0pt5_L2M3(double(L1L2_05_U1>0.5).*double(L1L2_05_U2<=0.5)==1 | double(L1L2_05_U1<=0.5).*double(L1L2_05_U2>0.5)==1)=0.1;

L1_0pt5_L2M = zeros(n,n,3);
L1_0pt5_L2M(:,:,1)=L1_0pt5_L2M1;
L1_0pt5_L2M(:,:,2) = L1_0pt5_L2M2;
L1_0pt5_L2M(:,:,3) = L1_0pt5_L2M3;

L1_0pt5_L2M = rescale_color_image(L1_0pt5_L2M);

%anisotropic
L1_M1=ones(n,n);
L1_M3=zeros(n,n);
L1_M1(double(ani_U1>0.5).*double(ani_U2>0.5)==1 | double(ani_U1>0.5).*double(ani_U2<=0.5)==1)=0.7;
L1_M2 = (double(ani_U1>0.5).*double(ani_U2>0.5)==1)*0.1+ (double(ani_U1<=0.5).*double(ani_U2>0.5)==1)*0.7;
L1_M3(double(ani_U1>0.5).*double(ani_U2<=0.5)==1 | double(ani_U1<=0.5).*double(ani_U2>0.5)==1)=0.1;

L1_M = zeros(n,n,3);
L1_M(:,:,1)=L1_M1;
L1_M(:,:,2) = L1_M2;
L1_M(:,:,3) = L1_M3;

L1_M = rescale_color_image(L1_M);

%isotropic
iso_M1=ones(n,n);
iso_M3=zeros(n,n);
iso_M1(double(iso_U1>0.5).*double(iso_U2>0.5)==1 | double(iso_U1>0.5).*double(iso_U2<=0.5)==1)=0.7;
iso_M2 = (double(iso_U1>0.5).*double(iso_U2>0.5)==1)*0.1+ (double(iso_U1<=0.5).*double(iso_U2>0.5)==1)*0.7;
iso_M3(double(iso_U1>0.5).*double(iso_U2<=0.5)==1 | double(iso_U1<=0.5).*double(iso_U2>0.5)==1)=0.1;

iso_M = zeros(n,n,3);
iso_M(:,:,1)=iso_M1;
iso_M(:,:,2) = iso_M2;
iso_M(:,:,3) = iso_M3;

iso_M = rescale_color_image(iso_M);

%compute ssim
ssim(L1_L2M, fg1)
ssim(L1_0pt5_L2M, fg1)
ssim(L1_M, fg1)
ssim(iso_M, fg1)

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