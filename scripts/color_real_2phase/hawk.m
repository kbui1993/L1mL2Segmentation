%read hawk image
fui8 = imread('hawk.jpg');

%convert to double
f = double(fui8);

%get size of image
[N,M,~] = size(f);

%rescale image
fg = rescale_color_image(f);


%set parameters
pm.outer_iter = 20;
pm.alpha = 1.0;
pm.lambda = 1.5;
pm.c = 1e-8;
pm.inner_iter = 300;
pm.tau = 1/6;
pm.sigma = 1/6;
pm.method = 'PDHG';


%set segmentation initialization
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


%plot segmentation
figure;
subplot(2,3,1); imagesc(fg); axis off; axis square; title('Original');
subplot(2,3,2); imagesc(fg); hold on; contour(double(L1_L2_u1>0.5), 'g'); axis off; axis square; title('L1-L2');
subplot(2,3,3); imagesc(fg); hold on; contour(double(L1_0pt5_L2_u1>0.5), 'g'); axis off; axis square; title('L1-0.5L2');
subplot(2,3,5); imagesc(fg); hold on; contour(double(L1_u1>0.5), 'g'); axis off; axis square; title('Anisotropic');
subplot(2,3,6); imagesc(fg); hold on; contour(double(iso_u1>0.5), 'g'); axis off; axis square; title('Isotropic');
