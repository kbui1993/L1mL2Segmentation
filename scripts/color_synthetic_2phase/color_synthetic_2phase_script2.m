%This script runs experiments to compare the proposed four-phase method 
%against isotropic and anisotropic Chan Vese segmentation on synthetic
%RGB image.

%read image
fui8 = imread('shape.png');
fui8 = rgb2gray(fui8);

f = double(fui8);
[N,M] = size(f);

%create color version
f_color = zeros(385,385,3);
f_color(:,:,1) = double(f>0)*0.5;
f_color(:,:,2) = double(f>0)*0.9;
f_color(:,:,3) = double(f>0)*0.25;

%rescale image
fg = rescale_color_image(f_color);

fg = double(fg);

%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 10;
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
L1_L2M1=zeros(N,M);
L1_L2M2=zeros(N,M);
L1_L2M3=zeros(N,M);
L1_L2M1(double(L1_L2_u1>=0.5)==1)=0.5;
L1_L2M2(double(L1_L2_u1>=0.5)==1)=0.9;
L1_L2M3(double(L1_L2_u1>=0.5)==1)=0.25;


L1_L2M = zeros(N,M,3);
L1_L2M(:,:,1)=L1_L2M1;
L1_L2M(:,:,2) = L1_L2M2;
L1_L2M(:,:,3) = L1_L2M3;

%L1-0.5L2 reconstruction
L1_0pt5_L2M1=zeros(N,M);
L1_0pt5_L2M2=zeros(N,M);
L1_0pt5_L2M3=zeros(N,M);
L1_0pt5_L2M1(double(L1_0pt5_L2_u1>=0.5)==1)=0.5;
L1_0pt5_L2M2(double(L1_0pt5_L2_u1>=0.5)==1)=0.9;
L1_0pt5_L2M3(double(L1_0pt5_L2_u1>=0.5)==1)=0.25;

L1_0pt5_L2M = zeros(N,M,3);
L1_0pt5_L2M(:,:,1)=L1_0pt5_L2M1;
L1_0pt5_L2M(:,:,2) = L1_0pt5_L2M2;
L1_0pt5_L2M(:,:,3) = L1_0pt5_L2M3;

%anisotropic reconstruction
L1_M1=zeros(N,M);
L1_M2=zeros(N,M);
L1_M3=zeros(N,M);
L1_M1(double(L1_u1>=0.5)==1)=0.5;
L1_M2(double(L1_u1>=0.5)==1)=0.9;
L1_M3(double(L1_u1>=0.5)==1)=0.25;

L1_M = zeros(N,M,3);
L1_M(:,:,1)=L1_M1;
L1_M(:,:,2) = L1_M2;
L1_M(:,:,3) = L1_M3;

%isotropic reconstruction
iso_M1=zeros(N,M);
iso_M2=zeros(N,M);
iso_M3=zeros(N,M);
iso_M1(double(iso_u1>=0.5)==1)=0.5;
iso_M2(double(iso_u1>=0.5)==1)=0.9;
iso_M3(double(iso_u1>=0.5)==1)=0.25;

iso_M = zeros(N,M,3);
iso_M(:,:,1)=iso_M1;
iso_M(:,:,2) = iso_M2;
iso_M(:,:,3) = iso_M3;

%compute DICE
dice(sum(rescale_color_image(L1_L2M),3)/3, sum(fg,3)/3)
dice(sum(rescale_color_image(L1_0pt5_L2M),3)/3, sum(fg,3)/3)
dice(sum(rescale_color_image(L1_M),3)/3, sum(fg,3)/3)
dice(sum(rescale_color_image(iso_M),3)/3, sum(fg,3)/3)

%plot segmentation
figure;
subplot(2,3,1); imagesc(f_color); axis off; axis square; title('Original');
subplot(2,3,2); imagesc(f_color); hold on; contour(double(L1_L2_u1>0.5), 'b'); axis off; axis square; title('L1-L2');
subplot(2,3,3); imagesc(f_color); hold on; contour(double(L1_0pt5_L2_u1>0.5), 'b'); axis off; axis square; title('L1-0.5L2');
subplot(2,3,5); imagesc(f_color); hold on; contour(double(L1_u1>0.5), 'b'); axis off; axis square; title('Anisotropic');
subplot(2,3,6); imagesc(f_color); hold on; contour(double(iso_u1>0.5), 'b'); axis off; axis square; title('Isotropic');

