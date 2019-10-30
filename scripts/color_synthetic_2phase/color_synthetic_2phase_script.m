%This script runs experiments to compare the proposed two-phase method 
%against isotropic and anisotropic Chan Vese segmentation on synthetic
%RGB image.
color_synthetic_image;
f = M;
[N,M,~] = size(f);
fg = rescale_color_image(f);

%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 1;
pm.c = 0;
pm.inner_iter = 300;
pm.tau = 1/8;
pm.sigma = 1/8;
pm.method = 'PDHG';

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
%L1-L2 reconstruction
L1_L2M1=ones(n,n);
L1_L2M2=zeros(n,n);
L1_L2M3=zeros(n,n);
L1_L2M1(double(L1_L2_u1<0.5)==1)=0;
L1_L2M2(double(L1_L2_u1<0.5)==1)=90;
L1_L2M3(double(L1_L2_u1<0.5)==1)=150;


L1_L2M = zeros(n,n,3);
L1_L2M(:,:,1)=L1_L2M1;
L1_L2M(:,:,2) = L1_L2M2;
L1_L2M(:,:,3) = L1_L2M3;

%L1-0.5L2 reconstruction
L1_0pt5_L2M1=ones(n,n);
L1_0pt5_L2M2=zeros(n,n);
L1_0pt5_L2M3=zeros(n,n);
L1_0pt5_L2M1(double(L1_0pt5_L2_u1<0.5)==1)=0;
L1_0pt5_L2M2(double(L1_0pt5_L2_u1<0.5)==1)=90;
L1_0pt5_L2M3(double(L1_0pt5_L2_u1<0.5)==1)=150;

L1_0pt5_L2M = zeros(n,n,3);
L1_0pt5_L2M(:,:,1)=L1_0pt5_L2M1;
L1_0pt5_L2M(:,:,2) = L1_0pt5_L2M2;
L1_0pt5_L2M(:,:,3) = L1_0pt5_L2M3;

%anisotropic reconstruction
L1_M1=ones(n,n);
L1_M2=zeros(n,n);
L1_M3=zeros(n,n);
L1_M1(double(L1_u1<0.5)==1)=0;
L1_M2(double(L1_u1<0.5)==1)=90;
L1_M3(double(L1_u1<0.5)==1)=150;

L1_M = zeros(n,n,3);
L1_M(:,:,1)=L1_M1;
L1_M(:,:,2) = L1_M2;
L1_M(:,:,3) = L1_M3;

%isotropic reconstruction
iso_M1=ones(n,n);
iso_M2=zeros(n,n);
iso_M3=zeros(n,n);
iso_M1(double(iso_u1<0.5)==1)=0;
iso_M2(double(iso_u1<0.5)==1)=90;
iso_M3(double(iso_u1<0.5)==1)=150;

iso_M = zeros(n,n,3);
iso_M(:,:,1)=iso_M1;
iso_M(:,:,2) = iso_M2;
iso_M(:,:,3) = iso_M3;

%compute ssim
ssim(L1_L2M, f)
ssim(L1_0pt5_L2M, f)
ssim(L1_M, f)
ssim(iso_M, f)

%plot segmentation
figure;
subplot(2,3,1); imagesc(f); axis off; axis square; title('Original');
subplot(2,3,2); imagesc(f); hold on; contour(double(L1_L2_u1>0.5), 'k'); axis off; axis square; title('L1-L2');
subplot(2,3,3); imagesc(f); hold on; contour(double(L1_0pt5_L2_u1>0.5), 'k'); axis off; axis square; title('L1-0.5L2');
subplot(2,3,5); imagesc(f); hold on; contour(double(L1_u1>0.5), 'k'); axis off; axis square; title('Anisotropic');
subplot(2,3,6); imagesc(f); hold on; contour(double(iso_u1>0.5), 'k'); axis off; axis square; title('Isotropic');
